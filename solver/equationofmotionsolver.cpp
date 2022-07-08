#include "equationofmotionsolver.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "spdlog/spdlog.h"
#include <Eigen/Core>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl.h"



void EquationOfMotionSolver::ClearAndResize(std::size_t N_)
{
    this->N=N_;
    if(cval.size() < N*dofs)
    {
        cval.resize(N*dofs*2);
        sln.resize(N*dofs*2);
    }

    if(csr_rows.size() < N+1) csr_rows.resize(N*2+1);

    std::fill(cval.begin(), cval.begin()+N*dofs, 0);

    if(rows_neighbors.size()<N)
        rows_neighbors.resize(N*2);

    for(unsigned i=0;i<N;i++) rows_neighbors[i].clear();
}

void EquationOfMotionSolver::AddNNZEntry(int row, int column)
{
    if(row < 0 || column < 0) return; // the element does not belong in the matrix
    if(row > column) std::swap(row,column);    // enforce upper-triangular matrix
    if((unsigned)row >= N) {
        spdlog::critical("trying to insert an element beyond the matrix size; row {}; N {}",row, N);
        throw std::runtime_error("trying to insert an element beyond the matrix size");
    }
    rows_neighbors[row].push_back(column);
}

void EquationOfMotionSolver::AddEntriesToStructure(const int* idx_begin, const int* idx_end)
{
    for(auto iter=(idx_begin+1); iter!=idx_end; ++iter)
        for(auto j=idx_begin; j!=iter; ++j)
            AddNNZEntry(*iter,*j);
}


void EquationOfMotionSolver::CreateStructure()
{
    // CREATE STRUCTURE ARRAYS

    // sort the neighbor list of each row
    int nnz_count = 0;
#pragma omp parallel for reduction(+:nnz_count)
    for(int i=0;i<(int)N;i++)
    {
        auto &rn = rows_neighbors[i];
        rn.push_back(i);    // add diagonal entry
        std::sort(rn.begin(),rn.end());
        rn.resize(std::distance(rn.begin(),std::unique(rn.begin(), rn.end())));     // remove duplicates
        nnz_count += rn.size();
    }

    csr_rows[N] = nnz = nnz_count;
    if(csr_cols.size() < nnz) csr_cols.resize(nnz*2);
    nnz*=dofssq;

    // ensure that the arrays are of sufficient sizes
    if(qoval.size() < nnz) qoval.resize(nnz*2);
    std::fill(qoval.begin(), qoval.begin()+nnz, 0);

    // enumerate entries
    unsigned count=0;
    for(unsigned row=0;row<N;row++)
    {
        csr_rows[row] = count;

        auto &rn = rows_neighbors[row];

        for(unsigned int const &local_column : rn)
        {
            if(row > local_column)
            {
                spdlog::critical("matrix is not upper-triangular; row {}; column {}",row, local_column);
                throw std::runtime_error("matrix is not upper-triangular");
            }
            csr_cols[count] = local_column;
            count++;
        }
    }

    if(nnz != dofssq*count)
    {
        spdlog::critical("csr_rows[{}]=={}, whereas count=={}",N,nnz,count);
        throw std::runtime_error("nnz != count");
    }
}

unsigned EquationOfMotionSolver::get_offset(const int row, const int column) const
{

    int col_offset_begin = csr_rows[row];
    int col_offset_end = csr_rows[row+1];

    const int *start_pt = &csr_cols[col_offset_begin];
    const int *end_pt = &csr_cols[col_offset_end];

    auto it = std::lower_bound(start_pt,end_pt,column);
    if(it == end_pt || *it != column)
    {
        spdlog::critical("offset for ({},{}) not found", row, column);
        throw std::runtime_error("EquationOfMotionSolver::get_offset(): (i,j) index not found");
    }
    unsigned offset = std::distance(start_pt,it)+col_offset_begin;
    offset*=dofssq;

    return offset;
}

void EquationOfMotionSolver::AddToQ(const int row, const int column, const double *v)
{
    if (row < 0 || column < 0 || row > column) return;
    else if((unsigned)row >= N || (unsigned)column >= N) throw std::runtime_error("AddToQ: out of range");
    int offset = get_offset(row,column);

    for(unsigned i=0;i<dofssq;i++)
    {
#pragma omp atomic
        qoval[offset+i] += v[i];
    }
}

void EquationOfMotionSolver::AddToC(const int idx, const double *v)
{
    if(idx < 0) return;
    if((unsigned)idx >= N) throw std::runtime_error("AddToC: index out of range");

    for(unsigned i=0;i<dofs;i++)
    {
#pragma omp atomic
        cval[idx*dofs+i] += v[i];
    }
}

void EquationOfMotionSolver::AddToEquation(const double *lE, const double *qE, const std::initializer_list<int> ids)
{
    unsigned n = ids.size();

    Eigen::Map<const Eigen::MatrixXd> levec(lE, 1, n*dofs);
    Eigen::Map<const Eigen::MatrixXd> qevec(qE, n*dofs, n*dofs);

    unsigned idx_i=0;
    for(auto iter_row=ids.begin(); iter_row!=ids.end(); ++iter_row,++idx_i)
    {
        int row = *iter_row;
        if(row < 0) continue;

        Eigen::Matrix<double, 1, dofs> vec = levec.block(idx_i*dofs,0, dofs,1);
        AddToC(row, vec.data());

        unsigned idx_j=0;
        for(auto iter_col=ids.begin(); iter_col!=ids.end(); ++iter_col,++idx_j)
        {
            int col = *iter_col;
            if(col < 0) continue;
            Eigen::Matrix<double,dofs,dofs> mat = qevec.block(idx_i*dofs, idx_j*dofs, dofs, dofs);
            mat.transposeInPlace();
            AddToQ(row,col,mat.data());
        }
    }
}


bool EquationOfMotionSolver::Solve()
{
    //SaveLinearSystem();   // for testing

    int n = N;
    MKL_INT mtype = -2;       // Real symmetric matrix: -2;  real unsymmetric: 11
    MKL_INT nrhs = 1;     // Number of right hand sides.
    void *pt[64] = {};
    MKL_INT iparm[64] = {};
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT idum;
    iparm[0] = 1;       // No solver default
    iparm[1] = 3;       // Fill-in reordering from METIS (was 2)
    iparm[3] = 0;       // No iterative-direct algorithm
    iparm[4] = 0;       // No user fill-in reducing permutation
    iparm[5] = 0;       // Write solution into x
    iparm[6] = 0;       // Not in use
    iparm[7] = 0;       // Max numbers of iterative refinement steps
    iparm[8] = 0;
    iparm[9] = 8;       // Perturb the pivot elements with 1E-iparm[9];
    iparm[10] = 0;      // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 0;
    iparm[12] = 0;      // Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy
    iparm[13] = 0;      // Output: Number of perturbed pivots
    iparm[14] = 0;
    iparm[15] = 0;
    iparm[16] = 0;
    iparm[17] = 1;      // 1 - disable report; Output: Number of nonzeros in the factor LU
    iparm[18] = 1;		// 1- disable report; output number of operations
    iparm[19] = 0;
    iparm[26] = 1;      // check matrix structure for errors; 0 - do not check
    iparm[27] = 0;      // 0 double; 1 single
    iparm[34] = 1;      // zero-base index
    iparm[36] = dofs;    // BSR with block size DOFS
    maxfct = 1;
    mnum = 1;
    msglvl = 0; // use 1 for verbose output
    error = 0;
    phase = 13;

    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, qoval.data(), csr_rows.data(), csr_cols.data(),
            &idum, &nrhs, iparm, &msglvl, cval.data(), sln.data(), &error);

    if(error != 0) throw std::runtime_error("MKL solver error");

    phase = -1; //clean up
    double ddum;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, csr_rows.data(), csr_cols.data(),
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if(error != 0) throw std::runtime_error("MKL solver error");

    double solution_norm_local = 0;
#pragma omp parallel for reduction(+:solution_norm_local)
    for(int i=0;i<(int)N*dofs;i++) solution_norm_local += (double)(sln[i]*sln[i]);

    solution_norm = sqrt(solution_norm_local);

    solveCount++;
    return true;
}


void EquationOfMotionSolver::GetTentativeResult(int idx, Eigen::Vector3d &vec)
{
    for(int i=0;i<dofs;i++) vec[i] = sln[idx*dofs+i];
}

void EquationOfMotionSolver::SaveLinearSystem()
{
    std::ostringstream fileName;

    fileName << "matrix_" << std::setfill('0') << std::setw(5) << solveCount << ".h";

    std::ofstream myfile;
    myfile.open(fileName.str());

    myfile << std::setprecision(9) << std::scientific;
    myfile << "constexpr int n = " << N << ";\n";
    myfile << "constexpr int nnz = " << nnz << ";\n";
    myfile << "double cval[] = {";
    for(unsigned i=0;i<N*dofs;i++) {
        myfile << cval[i];
        if(i!=N*dofs-1) myfile << ",";
        if(i%10==0) myfile << "\n";
    }
    myfile << "};\n";

    myfile << "double qval[] = {";
    for(unsigned i=0;i<nnz;i++) {
        myfile << qoval[i];
        if(i<nnz-1) myfile << ",";
        if(i%10==0) myfile << "\n";
    }
    myfile << "};\n";

    myfile << "int csr_rows[] = {";
    for(unsigned i=0;i<=N;i++) {
        myfile << csr_rows[i];
        if(i!=N) myfile << ",";
        if(i%30==0) myfile << "\n";
    }
    myfile << "};\n";

    myfile << "int csr_cols[] = {";
    for(unsigned i=0;i<nnz/dofssq;i++) {
        myfile << csr_cols[i];
        if(i<nnz/dofssq-1) myfile << ",";
        if(i%30==0) myfile << "\n";
    }
    myfile << "};\n";
    myfile.close();
}

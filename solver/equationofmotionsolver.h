#ifndef EQUATIONOFMOTIONSOLVER_H
#define EQUATIONOFMOTIONSOLVER_H

#include <tbb/concurrent_vector.h>

#include <map>
#include <unordered_map>
#include <initializer_list>
#include <vector>
#include <memory>

#include <Eigen/Core>


class EquationOfMotionSolver
{
public:

    void ClearAndResize(std::size_t N);     // size N must be set; return execution time

    void AddEntriesToStructure(const int* idx_begin, const int* idx_end); // insert nxn matrix of indices of non-zero entries
    void CreateStructure();

    // add values to non-zero elements
    void AddToEquation(const double *linearEntries, const double *quadraticEntries, const std::initializer_list<int> ids);

    constexpr static unsigned dofs = 5; // number of degrees of freedom per node (size of per-node blocks)
    constexpr static unsigned dofssq = dofs*dofs;

    // solver
    double solution_norm;
    bool Solve();   // true if successful
    void GetTentativeResult(int idx, Eigen::Vector3d &vec);  // solution => convenient vector form

private:
    std::vector<double> qoval;
    std::vector<int> csr_rows, csr_cols;

    // linear term
    std::vector<double> cval, sln;

    unsigned N;      // number of variables (divided by DOF)
    unsigned nnz;    // number of non-zero entries in Q (lower triangle)

    std::vector<tbb::concurrent_vector<unsigned>> rows_neighbors;  // list of indices of nz-columns per row

    void AddNNZEntry(int row, int column);    // reserve non-zero positions one-by-one (thread-safe)

    void ResizeRows();

    void AddToQ(const int row, const int column, const double *v);
    void AddToC(const int idx, const double *v);
    unsigned get_offset(const int row, const int column) const;

    int solveCount = 0;
    void SaveLinearSystem();    // for testing
};

#endif // EQUATIONOFMOTIONSOLVER_H

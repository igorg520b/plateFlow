#include "model.h"
#include "node.h"
#include "element.h"
#include "cohesivezone.h"

#include <QDebug>
#include <QThread>

#include <spdlog/spdlog.h>

icy::Model::Model() { Reset(); }

icy::Model::~Model() { HDF5Close(); }

void icy::Model::Reset()
{
    HDF5Close();

    currentStepDetails.Reset();
    replayMode = false;
    abortRequested = false;

    stepHistory.clear();
    stepHistory.reserve(2000);
    stepHistory.push_back(currentStepDetails);
    mesh.Reset();
    prms.Reset();
}

void icy::Model::Prepare()
{
    spdlog::info("icy::Model::Prepare()");
    abortRequested = false;
    if(!replayMode) stepHistory.resize(currentStepDetails.stepNumber+1);    // trim the list of saved steps
}

bool icy::Model::Step()
{
    if(replayMode)
    {
        std::unique_lock<std::mutex> mlock(replayMutex);
        cv_.wait(mlock);
        spdlog::info("replay time {}", replayTime);
        GoToSpecificTime(replayTime);
        replayTime += prms.p.VideoTimeStep;
        if(replayTime > stepHistory.back().time)
        {
            spdlog::info("replay done");
            GoToStep(0);
            replayMode = false;
            return false;
        }
        return true;
    }
    else
    return SimulationSingleStep();
}

void icy::Model::RequestAbort() { abortRequested = true; }

void icy::Model::Import(std::string fileName)
{
    Reset();
    mesh.Import(fileName);
}

void icy::Model::GoToStep(unsigned step)
{
    spdlog::info("icy::Model::GoToStep {}",step);
    if(step == 0 && file == nullptr)
    {
        currentStepDetails.Reset();
        mesh.GoToStep0();
    }
    else
    {
        if(file == nullptr) return;
        mesh.HDF5LoadStep(file, step); // load mesh state from HDF5
        currentStepDetails = stepHistory[step];
    }
    PositionKinematicObjects();
    ComputeVisualizedValues();
}

void icy::Model::Trim()
{
    stepHistory.resize(currentStepDetails.stepNumber+1);
    if(file != nullptr) icy::StepInfo::HDF5Trim(currentStepDetails.stepNumber, file);
}


void icy::Model::GoToSpecificTime(double whichTime)
{
    icy::StepInfo si_tmp;
    si_tmp.time = whichTime;
    auto it = std::lower_bound(stepHistory.begin(),stepHistory.end(),si_tmp,
                           [](const icy::StepInfo &s1, const icy::StepInfo &s2){return s1.time < s2.time;});
    unsigned offset = std::distance(stepHistory.begin(),it);
    double blendingParam;
    if(offset == stepHistory.size()-1) { offset--; blendingParam=1.;}
    else
    {
        double t1 = stepHistory[offset].time;
        double t2 = stepHistory[offset+2].time;
        blendingParam = (whichTime-t1)/(t2-t1);
    }
    mesh.HDF5LoadStepLERP(file, offset, blendingParam);
    spdlog::info("icy::Model::GoToSpecificTime {}; btw {} and {}; max {}", whichTime, offset, offset+1,stepHistory.size());
}







bool icy::Model::SimulationSingleStep()
{


    constexpr int colWidth = 12;

    int attempt = 0;
    bool converges=false;
    bool sln_okay;   // false if a CZ loads too fast or if solver cannot solve
    // gradually increase the time step
    constexpr double decayFactor = 1.1;
    double &tsf = currentStepDetails.tentativeStepFactor;
    tsf = std::clamp(tsf*decayFactor, 0., 1.);

    try
    {
        do
        {
            InitialGuess();
            spdlog::info("\n|{1:-^{0}}|{2:-^{0}}|{3:-^{0}}|{4:-^{0}}|{5:-^{0}}|{6:-^{0}}| ST {7:>}-{8:<2}",
                         colWidth, " it "," sln "," tsf "," ra ", " brdlst ", " coll ", currentStepDetails.stepNumber, attempt);

            double first_solution_norm = 0;
            int iter = 0;

            do
            {
                if(abortRequested) return false;
                sln_okay = AssembleAndSolve();

                double ratio = 0;
                if(iter == 0)  first_solution_norm = eqOfMotion.solution_norm;
                else if(first_solution_norm >= prms.p.ConvergenceCutoff) ratio = eqOfMotion.solution_norm/first_solution_norm;
                converges = (eqOfMotion.solution_norm < prms.p.ConvergenceCutoff ||
                             (ratio > 0 && ratio < prms.p.ConvergenceEpsilon));
                spdlog::info("|{1: ^{0}d}|{2: ^{0}.3e}|{3: ^{0}.3e}|{4: ^{0}.3e}|{5: ^{0}d}|{6: ^{0}d}|",
                             colWidth, iter, eqOfMotion.solution_norm, timeStepFactor_tentative,
                             ratio, mesh.broadlist.size(), mesh.pep_vec.size());
                iter++;
            } while(sln_okay && iter < prms.p.MaxIter && (iter < prms.p.MinIter || !converges));
            spdlog::info("|{0:-^{1}}|{0:-^{1}}|{0:-^{1}}|{0:-^{1}}|{0:-^{1}}|{0:-^{1}}|","",colWidth);

            if(!(sln_okay && converges))
            {
                spdlog::info("discarding attempt {}; sln_okay {}; converges {}",attempt,sln_okay,converges);
                attempt++;
                timeStepFactor_tentative*=0.5;
            }
            if(attempt > 30) throw std::runtime_error("Model::Step() could not solve");
            if(abortRequested) return false;

        } while (!(converges && sln_okay));

        AcceptTentativeValues();
    }
    catch(...)
    {
        HDF5Close();
        spdlog::critical("error occurred during Step()");
        return false;
    }

    return(currentStep < prms.p.MaxSteps && mesh.CountPinnedNodes() < mesh.nodes.size()/2);
}




bool icy::Model::AssembleAndSolve()
{
    unsigned nElems = mesh.elems.size();
    unsigned nNodes = mesh.nodes.size();
    unsigned nCZs = mesh.czs.size();

    // assign sequential indices to free nodes
    unsigned freeNodeCount = 0;
    for(Node *nd : mesh.nodes) nd->eqId = -1;

    for(unsigned i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.nodes[i];
        if(!nd->pinned) nd->eqId = freeNodeCount++;
        else nd->eqId = -1;
    }
    eqOfMotion.ClearAndResize(freeNodeCount);

#pragma omp parallel for
    for(int i=0;i<nElems;i++)
        mesh.elems[i]->AddToSparsityStructure(eqOfMotion);

    if(prms.p.EnableCZs)
    {
#pragma omp parallel for
        for(int i=0;i<nCZs;i++) mesh.czs[i]->AddToSparsityStructure(eqOfMotion);
    }

    if(prms.p.EnableCollisions)
    {
        mesh.UpdateTree();
        mesh.DetectContactPairs();
        mesh.AddToSparsityStructure(eqOfMotion);
    }


    eqOfMotion.CreateStructure();



    // ASSEMBLE

    bool mesh_inversion_detected = false;
    double h = timeStepFactor_tentative*prms.p.InitialTimeStep;
#pragma omp parallel for
    for(int i=0;i<nElems;i++)
    {
        if(mesh_inversion_detected) continue;
        bool elem_entry_ok = mesh.elems[i]->ComputeEquationEntries(eqOfMotion, prms.p, h);
        if(!elem_entry_ok) mesh_inversion_detected = true;
    }

    if(mesh_inversion_detected)
    {
        spdlog::info("icy::Model::AssembleAndSolve(): mesh inversion");
        return false; // mesh inversion
    }

    if(prms.p.EnableCZs)
    {
        bool cz_loading_too_fast = false;
#pragma omp parallel for
        for(int i=0;i<nCZs;i++)
        {
            if(cz_loading_too_fast) continue;
            bool cz_entry_ok = mesh.czs[i]->ComputeEquationEntries(eqOfMotion, prms.p);
            if(!cz_entry_ok) cz_loading_too_fast = true;
        }
        if(cz_loading_too_fast)
        {
            spdlog::info("icy::Model::AssembleAndSolve(): cz loading too fast");
            return false;
        }
    }

    if(prms.p.LinearIndenter)
    {
#pragma omp parallel for
    for(int i=0;i<nNodes;i++)
        mesh.nodes[i]->ComputeEquationEntries_Linear(eqOfMotion, prms.p, -indenter_position+1.0+prms.p.R_indenter);
    }
    else
    {
#pragma omp parallel for
    for(int i=0;i<nNodes;i++)
        mesh.nodes[i]->ComputeEquationEntries(eqOfMotion, prms.p, -indenter_position+1.0+prms.p.R_indenter);
    }

    if(prms.p.EnableCollisions)
    {
        if(prms.p.LinearCollisions)
            mesh.ComputeEquationEntries_Linear(eqOfMotion, prms.p);
        else
            mesh.ComputeEquationEntries_Log(eqOfMotion, prms.p);
    }

    // solve
    bool solve_result = eqOfMotion.Solve();

    if(!solve_result)
    {
        spdlog::info("eqOfMotion.Solve() returns {}",solve_result);
        return false;
    }
    if(std::isnan(eqOfMotion.solution_norm))
    {
        spdlog::critical("AssembleAndSovle: eqOfMotion.solution_norm == {}",eqOfMotion.solution_norm);
        throw std::runtime_error("solution_norm is NaN");
    }

// pull into Node::xt
#pragma omp parallel for
    for(int i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.nodes[i];
        Eigen::Vector3d delta_x;
        if(!nd->pinned)
        {
            eqOfMotion.GetTentativeResult(nd->eqId, delta_x);
            nd->xt+=delta_x;
        }
    }

    return true;
}

void icy::Model::PositionIndenter()
{
    // set the position of the indenter
    double h = timeStepFactor_tentative*prms.p.InitialTimeStep;
    indenter_position = (simulationTime_current+h)*prms.p.indentation_rate;
}


void icy::Model::InitialGuess()
{
    double h = timeStepFactor_tentative*prms.p.InitialTimeStep;
    PositionIndenter();
    std::size_t nNodes = mesh.nodes.size();
#pragma omp parallel for
    for(int i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.nodes[i];
        if(nd==nullptr) { spdlog::critical("InitialGuess: nd[{}]==nullptr",i); throw std::runtime_error("InitialGuess");}
        if(!nd->pinned) nd->xt = nd->xn + h*nd->vn;
    }
}

void icy::Model::AcceptTentativeValues()
{
    unsigned nNodes = mesh.nodes.size();
    unsigned nCZs = mesh.czs.size();
    double h = timeStepFactor_tentative*p.InitialTimeStep;

    mesh.mutexMeshUpdate.lock();

#pragma omp parallel for
    for(int i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.nodes[i];
        if(nd->pinned) continue;
        nd->vn = (nd->xt - nd->xn)/h;
        nd->xn = nd->xt;
    }

    mesh.mutexMeshUpdate.unlock();

    ComputeVisualizedValues();

    simulationTime_current += timeStepFactor_tentative*prms.p.InitialTimeStep;
    timeStepFactor_current = timeStepFactor_tentative;

    // calculate indenter force in z-dierection
    double z_force = 0;
    for(unsigned i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.nodes[i];
        if(nd->pinned) continue;
        z_force += nd->F[2];
    }

    //spdlog::info("zforce {}",z_force);

    currentStep++;

    StepInfo si;
    si.step = currentStep;
    si.time = simulationTime_current;
    si.timeStepFactor = timeStepFactor_current;
    si.czs_damaged = damaged_czs;
    si.czs_inactive = inactive_czs;
    si.indenter_force_z = -z_force;
    si.nCollisions = mesh.pep_vec.size();
    stepHistory.push_back(si);

    HDF5SaveStep();

    if(currentStep %20) mesh.DetectOutOfBounds();


}

void icy::Model::ComputeVisualizedValues()
{
    unsigned nElems = mesh.elems.size();
#pragma omp parallel for
    for(int i=0;i<nElems;i++) mesh.elems[i]->ComputeVisualizedVariables(prms.p);
}



void icy::Model::HDF5SaveNew(std::string fileName)
{
    HDF5Close();
    file = new H5::H5File(fileName, H5F_ACC_TRUNC);
    mesh.HDF5SaveInitialState(file);
    prms.HDF5SaveNew(file);
    icy::StepInfo::HDF5CreateDS(file);
    this->fileName = fileName;
    stepHistory.back().HDF5Save(file);
    mesh.HDF5SaveCurrentStep(file, currentStepDetails.stepNumber);
}

void icy::Model::HDF5Open(std::string fileName)
{
    Reset();    // saves the previous file and clears geometry
    file = new H5::H5File(fileName, H5F_ACC_RDWR);

    spdlog::info("icy::Model::HDF5Open reading params");
    prms.HDF5Load(file);
    spdlog::info("icy::Model::HDF5Open reading mesh");
    mesh.HDF5LoadInitialState(file);

    // read stepHistory
    icy::StepInfo::HDF5Read(file, stepHistory);
    this->fileName = fileName;

    spdlog::info("icy::Model::HDF5Open done {}",this->fileName);
}

void icy::Model::HDF5SaveStep()
{
    if(file == nullptr) return;
    stepHistory.back().HDF5Save(file);
    mesh.HDF5SaveCurrentStep(file, currentStepDetails.stepNumber);
}


void icy::Model::HDF5Close()
{
    if(file!=nullptr)
    {
        prms.HDF5SaveExisting(file);
        delete file;
        file = nullptr;
        fileName = "";
    }
}





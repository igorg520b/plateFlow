#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <string>
#include <vector>
#include <algorithm>
#include <mutex>
#include <condition_variable>

#include <H5Cpp.h>

#include "parameters_sim.h"
#include "modelcontrollerinterface.h"
#include "equationofmotionsolver.h"

#include "floemesh.h"
#include "modelstepinfo.h"


namespace icy { class Model; }

class icy::Model : public ModelControllerInterface
{
public:
    Model();
    ~Model();
    void Reset();

    Params p;
    icy::FloeMesh mesh;
    EquationOfMotionSolver solver;

    icy::StepInfo currentStep;
    std::vector<icy::StepInfo> stepHistory;

    void Prepare() override;        // invoked once, at simulation start
    bool Step() override;           // either invoked by Worker or via GUI
    void RequestAbort() override;   // invoked from GUI

    void Import(std::string fileName);
    void GoToStep(unsigned step);   // GUI instructs to load a given step from a file
    void Trim();                    // remove all history beyond current step

    // replay / save video
    std::string fileName;           // for GUI

private:
    bool abortRequested;

    bool SimulationSingleStep();
    void InitialGuess();
    bool AssembleAndSolve();        // return true if solved
    void AcceptTentativeValues();   // return true if plastic deformation occurred
    void ComputeVisualizedValues();
    void PositionKinematicObjects();

    // HDF5
public:
    void HDF5SaveNew(std::string fileName);
    void HDF5Open(std::string fileName);

private:
    H5::H5File *file = nullptr;
    void HDF5Close();
    void HDF5SaveStep();
};

#endif // MESHCOLLECTION_H

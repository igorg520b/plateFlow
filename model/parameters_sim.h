#ifndef P_SIM_H
#define P_SIM_H

#include <Eigen/Core>
#include <iostream>

#include <H5Cpp.h>

// variables related to the formulation of the model

namespace icy { struct Params; }


struct icy::Params
{
    bool EnableFracture;

    // integrator / simulation
    double InitialTimeStep;
    int MaxSteps, MinIter, MaxIter;
    double ConvergenceEpsilon, ConvergenceCutoff;

    // properties, material parameters
    double Gravity, WaterDensity, IceDensity, PoissonsRatio, YoungsModulus;
    double Thickness;

    // fracture
    double FractureTractionThreshold;
    double FractureWeakeningCoeff;
    double CutoffCoefficient; // used to reduce the computation cost of max_traction
    double FractureAngleThreshold, FractureAreaThreshold;
    double SubsteppingTimestepFactor;
    int SubstepIterations, SubstepRadius, FractureMaxSubsteps;


    // misc
    double VideoTimeStep;

    Eigen::Matrix3d elasticityMatrix;        // this has to be pre-computed whenever Y and nu change
    Eigen::Matrix2d D_mats;



    //void HDF5SaveNew(H5::H5File *file);
    //void HDF5SaveExisting(H5::H5File *file);
    //void HDF5Load(H5::H5File *file);



    void Reset()
    {
        EnableFracture = true;

        // integrator / simulation
        InitialTimeStep = 0.05;

        MaxSteps = 1000;
        MinIter = 3;
        MaxIter = 8;

        ConvergenceEpsilon = 1e-2;
        ConvergenceCutoff = 1e-8;


        Gravity = 9.81;
        WaterDensity = 997;
        IceDensity = 910;
        PoissonsRatio = 0.3;
        YoungsModulus = 3.7e9;


        Thickness = 0.1;

        SubstepIterations = 3;
        SubstepRadius = 10;
        SubsteppingTimestepFactor = 1e-3;

        FractureWeakeningCoeff = 0.75;

        CutoffCoefficient = 0.4;

        FractureAngleThreshold = 10;
        FractureAreaThreshold = 1e-4;
        FractureTractionThreshold = 1e5;
        FractureMaxSubsteps = 1000;

        RecomputeMatrices();

        VideoTimeStep = 1./300.;
    }

    void RecomputeMatrices()
    {
        elasticityMatrix.setZero();
        double k = YoungsModulus / (1.0 - PoissonsRatio*PoissonsRatio);
        elasticityMatrix(0,0) = elasticityMatrix(1,1) = k;
        elasticityMatrix(0,1) = elasticityMatrix(1,0) = k*PoissonsRatio;
        elasticityMatrix(2,2) = YoungsModulus/((1.0 + PoissonsRatio)*2.0);
        D_mats.setIdentity();
        D_mats *= ((5.0/6.0)*YoungsModulus/((1.0 + PoissonsRatio)*2.0));
    }
};



#endif

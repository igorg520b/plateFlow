#include "parameters_sim.h"
#include <spdlog/spdlog.h>

/*

void icy::SimParams::HDF5SaveNew(H5::H5File *file)
{
    H5::CompType t(sizeof(icy::PureParams));

    t.insertMember("EnableCollisions", offsetof(icy::PureParams, EnableCollisions), H5::PredType::NATIVE_HBOOL);
    t.insertMember("EnableCZs", offsetof(icy::PureParams, EnableCZs), H5::PredType::NATIVE_HBOOL);
    t.insertMember("LinearCollisions", offsetof(icy::PureParams, LinearCollisions), H5::PredType::NATIVE_HBOOL);
    t.insertMember("LinearIndenter", offsetof(icy::PureParams, LinearIndenter), H5::PredType::NATIVE_HBOOL);

    t.insertMember("YoungsModulus", offsetof(icy::PureParams, YoungsModulus), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("MaxSteps", offsetof(icy::PureParams, MaxSteps), H5::PredType::NATIVE_INT);
    t.insertMember("InitialTimeStep", offsetof(icy::PureParams, InitialTimeStep), H5::PredType::NATIVE_DOUBLE);

    t.insertMember("Gravity", offsetof(icy::PureParams, Gravity), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("Density", offsetof(icy::PureParams, Density), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("Kappa", offsetof(icy::PureParams, Kappa), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("PoissonsRatio", offsetof(icy::PureParams, PoissonsRatio), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("ConvergenceEpsilon", offsetof(icy::PureParams, ConvergenceEpsilon), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("ConvergenceCutoff", offsetof(icy::PureParams, ConvergenceCutoff), H5::PredType::NATIVE_DOUBLE);

    t.insertMember("MinIter", offsetof(icy::PureParams, MinIter), H5::PredType::NATIVE_INT);
    t.insertMember("MaxIter", offsetof(icy::PureParams, MaxIter), H5::PredType::NATIVE_INT);

    t.insertMember("R_indenter", offsetof(icy::PureParams, R_indenter), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("indentation_rate", offsetof(icy::PureParams, indentation_rate), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("ind_Kappa", offsetof(icy::PureParams, ind_Kappa), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("ind_dhat", offsetof(icy::PureParams, ind_dhat), H5::PredType::NATIVE_DOUBLE);

    t.insertMember("cz_alpha", offsetof(icy::PureParams, cz_alpha), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("cz_beta", offsetof(icy::PureParams, cz_beta), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("cz_lambda_n", offsetof(icy::PureParams, cz_lambda_n), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("cz_lambda_t", offsetof(icy::PureParams, cz_lambda_t), H5::PredType::NATIVE_DOUBLE);

    t.insertMember("cz_phi_n", offsetof(icy::PureParams, cz_phi_n), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("cz_phi_t", offsetof(icy::PureParams, cz_phi_t), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("cz_sigma_max", offsetof(icy::PureParams, cz_sigma_max), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("cz_tau_max", offsetof(icy::PureParams, cz_tau_max), H5::PredType::NATIVE_DOUBLE);

    t.insertMember("VideoTimeStep", offsetof(icy::PureParams, VideoTimeStep), H5::PredType::NATIVE_DOUBLE);

    hsize_t dim[1] = {1};
    H5::DataSpace dspace(1,dim);
    H5::DataSet dataset = file->createDataSet("/Initial/Params", t, dspace);
    dataset.write(&p, t);
}


void icy::SimParams::HDF5SaveExisting(H5::H5File *file)
{
    H5::DataSet dataset = file->openDataSet("/Initial/Params");
    H5::CompType t = dataset.getCompType();
    dataset.write(&p, t);
}

void icy::SimParams::HDF5Load(H5::H5File *file)
{
    H5::DataSet dataset = file->openDataSet("/Initial/Params");
    H5::CompType t = dataset.getCompType();

    dataset.read(&p, t);
    p.RecomputeCZParams();
    p.RecomputeLamdaAndMu();
    Q_EMIT propertyChanged();
    spdlog::info("icy::SimParams::HDF5Load done");
}

*/


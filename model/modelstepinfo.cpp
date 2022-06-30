#include "modelstepinfo.h"
#include <QDebug>

void icy::StepInfo::HDF5Save(H5::H5File *file)
{
    H5::DataSet ds = file->openDataSet("Steps");

    hsize_t newdims[1] = {(hsize_t)stepNumber+1};
    ds.extend(newdims);

    H5::DataSpace dsp = ds.getSpace();
    hsize_t offset[1] = {(hsize_t)stepNumber};
    hsize_t dims[1] = {1};
    dsp.selectHyperslab(H5S_SELECT_SET, dims, offset);

    H5::DataSpace memspace(1,dims);

    H5::CompType t = ds.getCompType();
    ds.write(this, t, memspace, dsp);
}

void icy::StepInfo::HDF5Trim(int count, H5::H5File *file)
{
    H5::DataSet ds = file->openDataSet("Steps");
    hsize_t newdims[1] = {(hsize_t)count+1};
    ds.extend(newdims);
}


void icy::StepInfo::HDF5Read(H5::H5File *file, std::vector<icy::StepInfo> &vec)
{
    H5::DataSet ds = file->openDataSet("Steps");
    hsize_t dims[1];
    ds.getSpace().getSimpleExtentDims(dims);
    vec.resize(dims[0]);

    H5::CompType t = ds.getCompType();
    ds.read(vec.data(), ds.getCompType());
}

void icy::StepInfo::HDF5CreateDS(H5::H5File *file)
{
    H5::CompType t(sizeof(icy::StepInfo));
    t.insertMember("stepNumber", offsetof(icy::StepInfo, stepNumber), H5::PredType::NATIVE_INT);
    t.insertMember("time", offsetof(icy::StepInfo, time), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("currentStepFactor", offsetof(icy::StepInfo, currentStepFactor), H5::PredType::NATIVE_DOUBLE);
    t.insertMember("tentativeStepFactor", offsetof(icy::StepInfo, tentativeStepFactor), H5::PredType::NATIVE_DOUBLE);

    hsize_t dim_initial[1] = {0};
    hsize_t dim_max[1] = {H5S_UNLIMITED};
    H5::DataSpace dspace(1,dim_initial, dim_max);

    H5::DSetCreatPropList cparms;
    hsize_t chunk_dims[1] = {100};
    cparms.setChunk(1, chunk_dims);
    cparms.setDeflate(7);

    H5::DataSet dataset = file->createDataSet("Steps", t, dspace, cparms);

    qDebug() << "creating dataset Nodes";
    // create the dataset for nodal coordinates
    hsize_t dim_initial_nodes[2] = {0,6};
    hsize_t dim_max_nodes[2] = {H5S_UNLIMITED,6};
    H5::DataSpace dspace_nodes(2,dim_initial_nodes, dim_max_nodes);

    H5::DSetCreatPropList cparms_nodes;
    hsize_t chunk_dims_nodes[2] = {100000,6};
    cparms_nodes.setChunk(2, chunk_dims_nodes);
    cparms_nodes.setDeflate(7);
    H5::DataSet dataset_nodes = file->createDataSet("Nodes", H5::PredType::NATIVE_DOUBLE, dspace_nodes, cparms_nodes);

    // dataset for cohesive zone info
    hsize_t dim_initial_czs[2] = {0,8};
    hsize_t dim_max_czs[2] = {H5S_UNLIMITED,8};
    H5::DataSpace dspace_czs(2, dim_initial_czs, dim_max_czs);

    H5::DSetCreatPropList cparms_czs;
    hsize_t chunk_dims_czs[2] = {100000,8};
    cparms_czs.setChunk(2, chunk_dims_czs);
    cparms_czs.setDeflate(7);
    H5::DataSet dataset_czs = file->createDataSet("CZs", H5::PredType::NATIVE_DOUBLE, dspace_czs, cparms_czs);
}


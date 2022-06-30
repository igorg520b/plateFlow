#include "mesh.h"
#include "element.h"
#include "cohesivezone.h"
#include "czinsertiontool.h"
#include "additional_math.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>

#include <Eigen/Core>

#include "spdlog/spdlog.h"


icy::ConcurrentPool<icy::Node> icy::Mesh::NodeFactory(reserveConst);
icy::ConcurrentPool<icy::Element> icy::Mesh::ElementFactory(reserveConst);
icy::ConcurrentPool<icy::CohesiveZone> icy::Mesh::CZFactory(reserveConst);
icy::ConcurrentPool<icy::BVHN> icy::Mesh::BVHNRootsFactory(reserveConst);


icy::Mesh::Mesh()
{
    broadlist.reserve(reserveConst*10);
    pep_vec.reserve(reserveConst);
}

icy::Mesh::~Mesh()
{
    Reset();
}

void icy::Mesh::Reset()
{
    tree_update_counter = 0;
    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    CZFactory.release(czs);
    BVHNRootsFactory.releaseAll();
    BVHN::BVHNFactory.releaseAll();
    pep_vec.clear();
    leaves.clear();
    grain_roots.clear();
    broadlist.clear();
    per_grain_bvhns.clear();
}

void icy::Mesh::GoToStep0()
{
    for(Node *nd : nodes)
    {
        nd->xn = nd->x0;
        nd->vn.setZero();
    }

    for(CohesiveZone *cz : czs)
    {
        for(int i=0;i<3;i++) cz->tmax[i]=cz->pmax[i]=0;
        cz->isActive = true;
        cz->isDamaged = false;
    }
}


icy::Node* icy::Mesh::AddNode()
{
    Node* nd = NodeFactory.take();
    nd->Reset();
    nd->globId = (int)nodes.size();
    nodes.push_back(nd);
    return nd;
}

icy::Element* icy::Mesh::AddElement()
{
    Element* elem = ElementFactory.take();
    elem->Reset();
    elem->elemId = (int)elems.size();
    elems.push_back(elem);
    return elem;
}

icy::CohesiveZone* icy::Mesh::AddCZ()
{
    CohesiveZone *cz = CZFactory.take();
    cz->Reset();
    czs.push_back(cz);
    return cz;
}


void icy::Mesh::LoadMSH(const std::string &fileName)
{
    Reset();
    gmsh::clear();
    gmsh::open(fileName);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, std::size_t> mtags; // gmsh nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(4, nodeTags, nodeCoords, parametricCoords);

    mutexMeshUpdate.lock();

    NodeFactory.release(nodes);
    ElementFactory.release(elems);
    nodes.clear();
    elems.clear();

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

        Node *nd = AddNode();
        mtags[tag] = nd->globId;
        nd->x0 = nd->xn = Eigen::Vector3d(nodeCoords[i*3+0], nodeCoords[i*3+1], nodeCoords[i*3+2]);
    }

    // GET ELEMENTS - per grain (entity)
    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,3);

    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> tetraTags, nodeTagsInTetra;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(4, tetraTags, nodeTagsInTetra,entityTag);

        for(std::size_t i=0;i<tetraTags.size();i++)
        {
            icy::Element *elem = AddElement();
            elem->grainId = j;
            for(int k=0;k<4;k++) elem->nds[k] = nodes[mtags.at(nodeTagsInTetra[i*4+k])];
        }
    }

    std::vector<std::size_t> lineTags, nodeTagsInLines;
    gmsh::model::mesh::getElementsByType(1, lineTags, nodeTagsInLines);
    std::cout << "number of lines: " << lineTags.size() << std::endl;


    // center mesh
    Eigen::Vector3d offset(0.5,0.5,0);
    spdlog::info("offset {}, {}",offset.x(),offset.y());

    for(Node *nd : nodes)
    {
        nd->x0 -= offset;
        nd->xn = nd->xt = nd->x0;
        nd->pinned = nd->x0.z() < 1e-7;
//        if(nd->x0.z() > 1.0 - 1e-7) nd->pinned = true; // for testing
    }


    CZInsertionTool czit;
    czit.InsertCZs(*this);

    for(icy::Element *elem : elems) elem->Precompute();     // Dm matrix and volume
    CreateLeaves();
    MarkIncidentFaces();

    mutexMeshUpdate.unlock();
    gmsh::clear();
    spdlog::info("icy::Mesh::LoadMSH completed");
}



void icy::Mesh::HDF5SaveInitialState(H5::H5File *file)
{
    file->createGroup("/Initial",55);

    // NODES
    hsize_t dims_nodes_initial[2] = {nodes.size(), 4};
    H5::DataSpace dataspace_nodes(2, dims_nodes_initial);

    H5::DSetCreatPropList cparms_nodes;
    hsize_t chunk_dims[2] = {3000, 4};
    cparms_nodes.setChunk(2, chunk_dims);
    cparms_nodes.setDeflate(9);
    H5::DataSet dataset_nodes = file->createDataSet("/Initial/Nodes", H5::PredType::NATIVE_DOUBLE, dataspace_nodes, cparms_nodes);

    std::vector<double> nds_buffer(4*nodes.size());
    for(std::size_t i=0;i<nodes.size();i++)
    {
        for(int j=0;j<3;j++)
            nds_buffer[i*4+j] = nodes[i]->x0[j];
        nds_buffer[i*4+3] = nodes[i]->pinned ? 1 : 0;
    }
    dataset_nodes.write(nds_buffer.data(), H5::PredType::NATIVE_DOUBLE);

    // ELEMENTS
    hsize_t dims_elems_initial[2] = {elems.size(), 5};
    H5::DataSpace dataspace_elems(2, dims_elems_initial);

    H5::DSetCreatPropList cparms_elems;
    hsize_t chunk_dims_elems[2] = {3000,5};
    cparms_elems.setChunk(2,chunk_dims_elems);
    cparms_elems.setDeflate(9);
    H5::DataSet dataset_elems = file->createDataSet("/Initial/Elems",H5::PredType::NATIVE_INT, dataspace_elems, cparms_elems);

    std::vector<int> elems_buffer(5*elems.size());
    for(std::size_t i=0;i<elems.size();i++)
    {
        for(int j=0;j<4;j++) elems_buffer[i*5+j] = elems[i]->nds[j]->globId;
        elems_buffer[i*5+4] = elems[i]->grainId;
    }
    dataset_elems.write(elems_buffer.data(), H5::PredType::NATIVE_INT);


    // COHESIVE ZONES
    hsize_t dims_czs_initial[2] = {czs.size(), 10};
    H5::DataSpace dataspace_czs(2, dims_czs_initial);

    H5::DSetCreatPropList cparms_czs;
    hsize_t chunk_dims_czs[2] = {3000,10};
    cparms_czs.setChunk(2,chunk_dims_czs);
    cparms_czs.setDeflate(9);
    H5::DataSet dataset_czs = file->createDataSet("/Initial/CZs",H5::PredType::NATIVE_INT,dataspace_czs,cparms_czs);

    std::vector<int> czs_buffer(10*czs.size());
    for(std::size_t i=0;i<czs.size();i++)
    {
        for(int j=0;j<6;j++) czs_buffer[i*10+j] = czs[i]->nds[j]->globId;
        czs_buffer[i*10+6] = czs[i]->elems[0]->elemId;
        czs_buffer[i*10+7] = czs[i]->elems[1]->elemId;
        czs_buffer[i*10+8] = czs[i]->faceIds[0];
        czs_buffer[i*10+9] = czs[i]->faceIds[1];
    }
    dataset_czs.write(czs_buffer.data(),H5::PredType::NATIVE_INT);
}

void icy::Mesh::HDF5LoadInitialState(H5::H5File *file)
{
    mutexMeshUpdate.lock();
    spdlog::info("icy::Mesh::HDF5LoadInitialState: reading nodes");

    // NODES
    H5::DataSet dataset_nodes = file->openDataSet("/Initial/Nodes");
    hsize_t dims_nodes[2];
    dataset_nodes.getSpace().getSimpleExtentDims(dims_nodes, NULL);

    int nNodes = dims_nodes[0];
    std::vector<double> buffer_nodes(dims_nodes[0]*dims_nodes[1]);
    dataset_nodes.read(buffer_nodes.data(),H5::PredType::NATIVE_DOUBLE);

    for(int i=0;i<nNodes;i++)
    {
        Node *nd = AddNode();
        nd->x0 = nd->xn = nd->xt = Eigen::Vector3d(buffer_nodes[i*4+0], buffer_nodes[i*4+1], buffer_nodes[i*4+2]);
        nd->pinned = buffer_nodes[i*4+3]==0 ? false : true;
    }


    spdlog::info("icy::Mesh::HDF5LoadInitialState: reading elems");
    // ELEMS
    H5::DataSet dataset_elems = file->openDataSet("/Initial/Elems");
    hsize_t dims_elems[2];
    dataset_elems.getSpace().getSimpleExtentDims(dims_elems, NULL);

    int nElems = dims_elems[0];
    std::vector<int> buffer_elems(dims_elems[0]*dims_elems[1]);
    dataset_elems.read(buffer_elems.data(),H5::PredType::NATIVE_INT);

    for(int i=0;i<nElems;i++)
    {
        icy::Element *elem = AddElement();
        for(int j=0;j<4;j++) elems[i]->nds[j] = nodes[buffer_elems[i*5+j]];
        elems[i]->grainId = buffer_elems[i*5+4];
    }


    spdlog::info("icy::Mesh::HDF5LoadInitialState: reading czs");
    // CZS
    H5::DataSet dataset_czs = file->openDataSet("/Initial/CZs");
    hsize_t dims_czs[2];
    dataset_czs.getSpace().getSimpleExtentDims(dims_czs, NULL);

    int nCZs = dims_czs[0];
    std::vector<int> buffer_czs(dims_czs[0]*dims_czs[1]);
    dataset_czs.read(buffer_czs.data(),H5::PredType::NATIVE_INT);

    for(int i=0;i<nCZs;i++)
    {
        icy::CohesiveZone *cz = AddCZ();
        for(int j=0;j<6;j++) cz->nds[j] = nodes[buffer_czs[i*10+j]];
        cz->elems[0] = elems[buffer_czs[i*10+6]];
        cz->elems[1] = elems[buffer_czs[i*10+7]];
        cz->faceIds[0] = buffer_czs[i*10+8];
        cz->faceIds[1] = buffer_czs[i*10+9];

        for(int k=0;k<3;k++)
            if(cz->nds[k]->x0 != cz->nds[k+3]->x0)
                throw std::runtime_error("cz insertion error");
    }

    for(icy::Element *elem : elems) elem->Precompute();     // Dm matrix and volume
    CreateLeaves();
    MarkIncidentFaces();

    mutexMeshUpdate.unlock();
}

void icy::Mesh::HDF5SaveCurrentStep(H5::H5File *file, unsigned step)
{
    // prepare buffers
    std::vector<double> nodes_buffer(nodes.size()*6);
    std::vector<double> czs_buffer(nodes.size()*8);

    mutexMeshUpdate.lock();
    for(unsigned i=0;i<nodes.size();i++)
        for(unsigned j=0;j<3;j++)
    {
        nodes_buffer[i*6+j] = nodes[i]->xn[j];
        nodes_buffer[i*6+3+j] = nodes[i]->vn[j];
    }

    for(unsigned i=0;i<czs.size();i++)
    {
        CohesiveZone *cz = czs[i];
        for(int j=0;j<3;j++)
        {
            czs_buffer[i*8+j] = cz->pmax[j];
            czs_buffer[i*8+3+j] = cz->tmax[j];
        }
        czs_buffer[i*8+6] = cz->isActive ? 1. : 0;
        czs_buffer[i*8+7] = cz->isDamaged ? 1. : 0;
    }
    mutexMeshUpdate.unlock();

    // write nodes
    H5::DataSet dsNodes = file->openDataSet("Nodes");

    hsize_t dims_nodes[2] = {(step+1)*nodes.size(),6};
    dsNodes.extend(dims_nodes);

    H5::DataSpace dsp = dsNodes.getSpace();
    hsize_t   offset_nodes[2] = {(step)*nodes.size(),0};
    hsize_t dims1[2] = {nodes.size(), 6};
    dsp.selectHyperslab(H5S_SELECT_SET, dims1, offset_nodes);

    H5::DataSpace memspace(2, dims1);
    dsNodes.write(nodes_buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, dsp);

    // write czs
    H5::DataSet dsCZs = file->openDataSet("CZs");
    hsize_t dims_czs[2] = {(step+1)*czs.size(),8};
    dsCZs.extend(dims_czs);

    H5::DataSpace dsp_czs = dsCZs.getSpace();
    hsize_t   offset_czs[2] = {(step)*czs.size(),0};
    hsize_t dims1_czs[2] = {czs.size(), 8};
    dsp_czs.selectHyperslab(H5S_SELECT_SET, dims1_czs, offset_czs);

    H5::DataSpace memspace_czs(2, dims1_czs);
    dsCZs.write(czs_buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace_czs, dsp_czs);
}

void icy::Mesh::HDF5LoadStep(H5::H5File *file, unsigned step)
{
    H5::DataSet dsNodes = file->openDataSet("Nodes");
    H5::DataSet dsCZs = file->openDataSet("CZs");

    // prepare buffers
    std::vector<double> nodes_buffer(nodes.size()*6);
    std::vector<double> czs_buffer(nodes.size()*8);

    // read nodes
    H5::DataSpace dsp = dsNodes.getSpace();
    hsize_t   offset_nodes[2] = {(step)*nodes.size(),0};
    hsize_t dims1[2] = {nodes.size(), 6};
    dsp.selectHyperslab(H5S_SELECT_SET, dims1, offset_nodes);

    H5::DataSpace memspace(2, dims1);
    dsNodes.read(nodes_buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, dsp);

    // read CZs
    H5::DataSpace dsp_czs = dsCZs.getSpace();
    hsize_t   offset_czs[2] = {(step)*czs.size(),0};
    hsize_t dims1_czs[2] = {czs.size(), 8};
    dsp_czs.selectHyperslab(H5S_SELECT_SET, dims1_czs, offset_czs);

    H5::DataSpace memspace_czs(2, dims1_czs);
    dsCZs.read(czs_buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace_czs, dsp_czs);

    // distribute buffers into nodes and csz
    mutexMeshUpdate.lock();
    for(std::size_t i=0;i<nodes.size();i++)
        for(unsigned j=0;j<3;j++)
        {
            nodes[i]->xn[j] = nodes_buffer[i*6+j];
            nodes[i]->vn[j] = nodes_buffer[i*6+3+j];
        }

    for(unsigned i=0;i<czs.size();i++)
    {
        CohesiveZone *cz = czs[i];
        for(int j=0;j<3;j++)
        {
            cz->pmax[j] = czs_buffer[i*8+j];
            cz->tmax[j] = czs_buffer[i*8+3+j];
        }
        cz->isActive = czs_buffer[i*8+6] == 0 ? false : true;
        cz->isDamaged = czs_buffer[i*8+7] == 0 ? false : true;
    }
    mutexMeshUpdate.unlock();
}


void icy::Mesh::HDF5LoadStepLERP(H5::H5File *file, unsigned step, double b)
{
    H5::DataSet dsNodes = file->openDataSet("Nodes");
    H5::DataSet dsCZs = file->openDataSet("CZs");

    // prepare buffers
    std::vector<double> nodes_buffer(nodes.size()*6*2);
    std::vector<double> czs_buffer(nodes.size()*8);

    // read nodes
    H5::DataSpace dsp = dsNodes.getSpace();
    hsize_t   offset_nodes[2] = {(step)*nodes.size(),0};
    hsize_t dims1[2] = {nodes.size()*2, 6};
    dsp.selectHyperslab(H5S_SELECT_SET, dims1, offset_nodes);

    H5::DataSpace memspace(2, dims1);
    dsNodes.read(nodes_buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, dsp);

    // read CZs
    H5::DataSpace dsp_czs = dsCZs.getSpace();
    hsize_t   offset_czs[2] = {(step)*czs.size(),0};
    hsize_t dims1_czs[2] = {czs.size(), 8};
    dsp_czs.selectHyperslab(H5S_SELECT_SET, dims1_czs, offset_czs);

    H5::DataSpace memspace_czs(2, dims1_czs);
    dsCZs.read(czs_buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace_czs, dsp_czs);

    // distribute buffers into nodes and csz
    mutexMeshUpdate.lock();
    for(std::size_t i=0;i<nodes.size();i++)
        for(unsigned j=0;j<3;j++)
        {
            nodes[i]->xn[j] = (1-b)*nodes_buffer[i*6+j] + b*nodes_buffer[(i+nodes.size())*6+j];
            nodes[i]->vn[j] = (1-b)*nodes_buffer[i*6+3+j] + b*nodes_buffer[(i+nodes.size())*6+3+j];
        }

    for(unsigned i=0;i<czs.size();i++)
    {
        CohesiveZone *cz = czs[i];
        for(int j=0;j<3;j++)
        {
            cz->pmax[j] = czs_buffer[i*8+j];
            cz->tmax[j] = czs_buffer[i*8+3+j];
        }
        cz->isActive = czs_buffer[i*8+6] == 0 ? false : true;
        cz->isDamaged = czs_buffer[i*8+7] == 0 ? false : true;
    }
    mutexMeshUpdate.unlock();
}

void icy::Mesh::MarkIncidentFaces()
{
    // initialize incident faces information in nodes and elements

    // (1) find exposed faces
    std::map<std::tuple<int,int,int>,icy::CZInsertionTool::Facet> facets;
    for(icy::Element *e : elems)
    {
        for(int k=0;k<4;k++)
        {
            std::tuple<int,int,int> key = icy::CZInsertionTool::Facet::make_key(
                        e->nds[Element::fi[k][0]],
                    e->nds[Element::fi[k][1]],
                    e->nds[Element::fi[k][2]]);
            icy::CZInsertionTool::Facet facet;
            facet.key = key;
            facet.elems[0] = e;
            facet.facet_idx[0] = k;
            auto result = facets.insert({key,facet});
            if(result.second == false)
            {
                icy::CZInsertionTool::Facet &f = result.first->second;
                f.elems[1] = e;
                f.facet_idx[1] = k;
            }
        }
    }

    // (2) distribute exposed faces to their incident nodes
    for(Node *nd : nodes) nd->incident_faces.clear();

    int count = 0;
    for(auto &kvp : facets)
    {
        icy::CZInsertionTool::Facet &f = kvp.second;
        if(f.elems[1] != nullptr) continue;
        Element *e = f.elems[0];
        unsigned k = f.facet_idx[0];

        uint32_t facet_code = (uint32_t)f.facet_idx[0] | (uint32_t)f.elems[0]->elemId << 2;
        for(int i=0;i<3;i++)
        {
            Node *nd = e->nds[Element::fi[k][i]];
            nd->incident_faces.push_back(facet_code);
        }
        count++;
    }

    // (3) per element - collect all exposed faces from element's 4 nodes
    for(Element *elem : elems)
    {
        auto &r = elem->elem_incident_faces;
        r.clear();
        for(int i=0;i<4;i++)
            for(uint32_t facet_code : elem->nds[i]->incident_faces)
                r.push_back(facet_code);
        std::sort(r.begin(),r.end());
        r.resize(std::distance(r.begin(),std::unique(r.begin(), r.end())));
    }
}



void icy::Mesh::CreateLeaves()
{
    tree_update_counter = 0;
    BVHNRootsFactory.release(leaves);
    BVHNRootsFactory.release(grain_roots);

    auto iter = std::max_element(elems.begin(),elems.end(),
                                 [](icy::Element *e1, icy::Element *e2) {return e1->grainId < e2->grainId;});
    int nGrains = (*iter)->grainId + 1;
    spdlog::info("nGrains {}",nGrains);

    unsigned avg_grain_elems = elems.size()/nGrains*1.5;    // number of elements in grains
    per_grain_bvhns.resize(nGrains);
    for(auto vec : per_grain_bvhns) { vec.clear(); vec.reserve(avg_grain_elems); }

    grain_roots.resize(nGrains);
    for(int i=0;i<nGrains;i++)
    {
        BVHN *grain_root = BVHNRootsFactory.take();
        grain_roots[i] = grain_root;
        grain_root->elem = nullptr;
        grain_root->isLeaf = false;
        grain_root->test_self_collision = false;
    }

    leaves.reserve(elems.size());
    for(Element *elem : elems)
    {
        BVHN *leaf_ccd = BVHNRootsFactory.take();
        leaves.push_back(leaf_ccd);
        leaf_ccd->elem = elem;
        leaf_ccd->isLeaf = true;
        per_grain_bvhns[elem->grainId].push_back(leaf_ccd);
    }
}



void icy::Mesh::UpdateTree()
{
    // update leafs
    int nLeaves = (int)leaves.size();
#pragma omp parallel for
    for(int i=0;i<nLeaves;i++)
        leaves[i]->Expand();

    if(tree_update_counter % 10 == 0)
    {
        BVHN::BVHNFactory.releaseAll();

        if(grain_roots.size()>1)
        {
            for(int i=0;i<grain_roots.size();i++) grain_roots[i]->Build(per_grain_bvhns[i],0);
            mesh_root.Build(grain_roots,0);
        }
        else
        {
            mesh_root.Build(per_grain_bvhns.front(),0);
        }
    }
    else
    {
        mesh_root.Update();
    }
    tree_update_counter++;
}


void icy::Mesh::DetectContactPairs()
{
    // BROAD PHASE of contact detection
    broadlist.clear();
    mesh_root.SelfCollide(broadlist);

    // NARROW PHASE - check if tetrahedra intersect (each node vs opposite tetrahedron)
    int nBroadList = (int)broadlist.size();

    tbb::concurrent_unordered_set<uint64_t> point_element_pairs;

#pragma omp parallel for
    for(int i=0;i<nBroadList;i++)
    {
        BVHN *bvhn1, *bvhn2;
        std::tie(bvhn1,bvhn2) = broadlist[i];
        Element *elem1 = bvhn1->elem;
        Element *elem2 = bvhn2->elem;
        Eigen::Matrix3d m1 = elem1->getSpecialInverseDs();
        Eigen::Matrix3d m2 = elem2->getSpecialInverseDs();
        for(int j=0;j<4;j++)
        {
            if(elem1->IsNodeInside(m1, elem2->nds[j]->xt))
                point_element_pairs.insert((uint64_t)elem1->elemId << 32 | (uint64_t)elem2->nds[j]->globId);
            if(elem2->IsNodeInside(m2, elem1->nds[j]->xt))
                point_element_pairs.insert((uint64_t)elem2->elemId << 32 | (uint64_t)elem1->nds[j]->globId);
        }
    }

    // for each intersection, find nearest exposed face (if no face, set 0xFFFFFFFF code)
    pep_vec.resize(point_element_pairs.size());
    std::copy(point_element_pairs.begin(),point_element_pairs.end(),pep_vec.begin());

    int nCollisions = (int)pep_vec.size();

#pragma omp parallel for
    for(int i=0;i<nCollisions;i++)
    {
        uint64_t code = pep_vec[i];
        unsigned node_idx = code & 0xFFFFFFFF;
        unsigned elem_idx = code >> 32;
        Element *elem = elems[elem_idx];
        Node *nd = nodes[node_idx];
        pep_vec[i] = (uint64_t)FindClosestFace(nd, elem) | (uint64_t)node_idx<<32;
    }
}


uint32_t icy::Mesh::FindClosestFace(Node *nd, Element *elem)
{
    if(elem->elem_incident_faces.size() == 0) return 0xFFFFFFFF;

    Eigen::Vector3d pt = nd->xt;
    uint32_t closest_face;
    double smallest_distance = DBL_MAX;

    for(uint32_t face_code : elem->elem_incident_faces)
    {
        Node *nds[3];
        FaceCode2NodeCoords(face_code, nds);
        Eigen::Vector3d v[3] = {nds[0]->xt,nds[1]->xt,nds[2]->xt};
        double dist = icy::AdditionalMath::dtn(pt, v);
        if(dist < smallest_distance)
        {
            smallest_distance = dist;
            closest_face = face_code;
        }
    }
    return closest_face;
}

void icy::Mesh::FaceCode2NodeCoords(uint32_t face_code, Node *nds[3]) const
{
    unsigned elem_idx = face_code >> 2;
    Element *fe = elems[elem_idx];
    unsigned k = face_code & 0x3;
    for(int i=0;i<3;i++) nds[i] = fe->nds[Element::fi[k][i]];
}



void icy::Mesh::AddToSparsityStructure(EquationOfMotionSolver &eq) const
{
    int nContacts = pep_vec.size();

#pragma omp parallel for
    for(int i=0;i<nContacts;i++)
    {
        uint64_t face_node_code = pep_vec[i];
        uint32_t node_idx = face_node_code >> 32;
        uint32_t face_code = face_node_code & 0xFFFFFFFF;
        if(face_code == 0xFFFFFFFF) continue; // nearest exposed face was not found (the intersecting element does not have one)

        Node *nds[3];
        FaceCode2NodeCoords(face_code, nds);

        Node *nd = nodes[node_idx];
        int idxs[] {nds[0]->eqId, nds[1]->eqId, nds[2]->eqId, nd->eqId};
        eq.AddEntriesToStructure(std::begin(idxs),std::end(idxs));
    }
}



void icy::Mesh::DetectOutOfBounds()
{
    std::unordered_set<int> grains_to_freeze;
    for(Element *elem : elems)
        if(!elem->nds[0]->pinned && elem->nds[0]->xn[2] < -1)
            grains_to_freeze.insert(elem->grainId);

    if(grains_to_freeze.size() == 0) return;

    for(Element *elem : elems)
    {
        if(grains_to_freeze.count(elem->grainId))
        {
            for(int i=0;i<4;i++)
            {
                Node *nd = elem->nds[i];
                nd->pinned = true;
                nd->xt = nd->xn;
                nd->vn.setZero();
            }
        }
    }
}

void icy::Mesh::ComputeEquationEntries_Log(EquationOfMotionSolver &eq, const PureParams &prms)
{
    constexpr double distanceEpsilonSqared = 1e-30;
    const double dhat = prms.dhat;

    int nContacts = pep_vec.size();

#pragma omp parallel for
    for(int i=0;i<nContacts;i++)
    {
        uint64_t face_node_code = pep_vec[i];
        uint32_t node_idx = face_node_code >> 32;
        uint32_t face_code = face_node_code & 0xFFFFFFFF;
        if(face_code == 0xFFFFFFFF) continue; // nearest exposed face was not found (the intersecting element does not have one)

        Node *nds[3];
        FaceCode2NodeCoords(face_code, nds);    // convert face_code to array of nodes
        Node *nd = nodes[node_idx];

        // calculate the distance squared, its gradient and Hessian
        Eigen::Matrix<double,12,1> g;
        Eigen::Matrix<double,12,12> H;
        Eigen::Vector3d pts[4] = {nd->xt, nds[0]->xt, nds[1]->xt, nds[2]->xt};
        double dsq = icy::AdditionalMath::pt_alt(pts, g, H);

        if (dsq < distanceEpsilonSqared) continue;

        // calculate grad/Hessian of distance
        double d = sqrt(dsq);
        if (d >= dhat) continue;
        Eigen::Matrix<double,12,1> gd = g/(2*d);
        Eigen::Matrix<double,12,12> Hd = H/(2*d) - g*g.transpose()/(4*d*dsq);

        // logarithmic potential
        double l = log((dhat-d)/dhat);
        //double b = -dsq*l;
        double bp = dsq/(dhat-d) -2*d*l;
        double bpp = dsq/((dhat-d)*(dhat-d)) + 4*d/(dhat-d) - 2*l;

        Eigen::Matrix<double,12,1> gb = gd*bp;
        Eigen::Matrix<double,12,12> Hb = bpp*gb*gb.transpose() + bp*Hd;

        // evaluate the forcing function
        Eigen::Matrix<double,12,1> rhs;
        Eigen::Matrix<double,12,12> lhs;

        double k = prms.Kappa/2;
        rhs = -k*gb;
        lhs = k*Hb;

        eq.AddToEquation(rhs.data(), lhs.data(), {nd->eqId, nds[0]->eqId, nds[1]->eqId, nds[2]->eqId});
    }
}


void icy::Mesh::ComputeEquationEntries_Linear(EquationOfMotionSolver &eq, const PureParams &prms)
{
    constexpr double distanceEpsilonSqared = 1e-20;
    int nContacts = pep_vec.size();

#pragma omp parallel for
    for(int i=0;i<nContacts;i++)
    {
        uint64_t face_node_code = pep_vec[i];
        uint32_t node_idx = face_node_code >> 32;
        uint32_t face_code = face_node_code & 0xFFFFFFFF;
        if(face_code == 0xFFFFFFFF) continue; // nearest exposed face was not found (the intersecting element does not have one)

        Node *nds[3];
        FaceCode2NodeCoords(face_code, nds);    // convert face_code to array of nodes
        Node *nd = nodes[node_idx];

        Eigen::Matrix<double,12,1> xt; // force, nodal coords
        xt << nd->xt, nds[0]->xt, nds[1]->xt, nds[2]->xt;

        // calculate the distance squared, its gradient and Hessian
        Eigen::Matrix<double,12,1> g;
        Eigen::Matrix<double,12,12> H;
        Eigen::Vector3d pts[4] = {nd->xt, nds[0]->xt, nds[1]->xt, nds[2]->xt};
        double dsq = icy::AdditionalMath::pt_alt(pts, g, H);

        if (dsq < distanceEpsilonSqared) continue;

        // evaluate the forcing function
        Eigen::Matrix<double,12,1> rhs;
        Eigen::Matrix<double,12,12> lhs;

        double k = prms.Kappa/2;
        rhs = -k*g;
        lhs = k*H;
        eq.AddToEquation(rhs.data(), lhs.data(), {nd->eqId, nds[0]->eqId, nds[1]->eqId, nds[2]->eqId});
    }
}



void icy::Mesh::ExportForAbaqus()
{
    std::ofstream s;
    s.open("ice_sample.py",std::ios_base::trunc|std::ios_base::out);
    s << std::setprecision(9);
    s << "from abaqus import *\n";
    s << "from abaqusConstants import *\n";
    s << "import mesh\n";
    s << "import regionToolset\n";
    s << "p = mdb.models['Model-1'].Part(name='MyPart1', dimensionality=THREE_D, type=DEFORMABLE_BODY)\n";

    s << "print \"importing nodes\" \n";

    for(Node *nd : nodes)
        s << "p.Node(coordinates=(" << nd->x0[0] << "," << nd->x0[1] << "," << nd->x0[2] << "))\n";

    s << "n = p.nodes\n";

    for(Element *e : elems)
        s << "p.Element(nodes=(n["<<e->nds[0]->globId<<"],n["<<e->nds[1]->globId<<
             "],n["<<e->nds[3]->globId<<"],n["<<e->nds[2]->globId<<"]), elemShape=TET4)\n";

    for(icy::CohesiveZone *c : czs)
        s << "p.Element(nodes=(n["<<c->nds[0]->globId<<
             "], n["<<c->nds[1]->globId<<
             "], n["<<c->nds[2]->globId<<
             "], n["<<c->nds[3]->globId<<
             "], n["<<c->nds[4]->globId<<
             "], n["<<c->nds[5]->globId<<"]), elemShape=WEDGE6)\n";

    s << "elemType_bulk = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)\n";
    s << "elemType_coh = mesh.ElemType(elemCode=COH3D6, elemLibrary=STANDARD)\n";

    // region1 - bulk elements

    s << "region1 = p.elements[0:"<<elems.size()<<"]\n";
    s << "p.setElementType(regions=(region1,), elemTypes=(elemType_bulk,))\n";
    s << "p.Set(elements=(region1,), name='Set-1-elems')\n";

    s << "region2cz = p.elements[" << elems.size() << ":" << elems.size()+czs.size() << "]\n";
    s << "p.setElementType(regions=(region2cz,), elemTypes=(elemType_coh,))\n";
    s << "p.Set(elements=(region2cz,), name='Set-2-czs')\n";

    // region2 - pinned nodes
    s << "region3pinned = (";
    for(Node *nd : nodes)
        if(nd->pinned)
            s << "p.nodes["<<nd->globId<<":"<<nd->globId+1<<"],";
    s << ")\n";
    s << "p.Set(nodes=region3pinned,name='Set3-pinned')\n";

    // create bulk material
    s << "mat1 = mdb.models['Model-1'].Material(name='Material-1-bulk')\n";
    s << "mat1.Density(table=((900.0, ), ))\n";
    s << "mat1.Elastic(table=((1000000000.0, 0.3), ))\n";

    // cz material
    s << "mat2 = mdb.models['Model-1'].Material(name='Material-2-czs')\n";
    s << "mat2.Density(table=((1.0, ), ))\n";
    s << "mat2.MaxsDamageInitiation(table=((100000.0, 50000.0, 50000.0), ))\n";
    s << "mat2.maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((30.0, ), ))\n";
    s << "mat2.Elastic(type=TRACTION, table=((2000000000.0, 1000000000.0, 1000000000.0), ))\n";

    // sections
    s << "mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1-bulk', "
        "material='Material-1-bulk', thickness=None)\n";
    s << "mdb.models['Model-1'].CohesiveSection(name='Section-2-czs', "
        "material='Material-2-czs', response=TRACTION_SEPARATION, "
        "outOfPlaneThickness=None)\n";

    // section assignments
    s << "region = p.sets['Set-1-elems']\n";
    s << "p.SectionAssignment(region=region, sectionName='Section-1-bulk', offset=0.0, "
        "offsetType=MIDDLE_SURFACE, offsetField='', "
        "thicknessAssignment=FROM_SECTION)\n";
    s << "region = p.sets['Set-2-czs']\n";
    s << "p = mdb.models['Model-1'].parts['MyPart1']\n";
    s << "p.SectionAssignment(region=region, sectionName='Section-2-czs', offset=0.0, "
        "offsetType=MIDDLE_SURFACE, offsetField='', "
        "thicknessAssignment=FROM_SECTION)\n";

    // indenter
    s << "s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)\n";
    s << "g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n";
    s << "s.setPrimaryObject(option=STANDALONE)\n";
    s << "s.ArcByCenterEnds(center=(0.0, 0.0), point1=(-0.05, 0.0), point2=(0.05, -0.0125), direction=COUNTERCLOCKWISE)\n";
    s << "p2 = mdb.models['Model-1'].Part(name='Part-2', dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)\n";
    //s << "p2 = mdb.models['Model-1'].parts['Part-2']\n";
    s << "p2.AnalyticRigidSurfExtrude(sketch=s, depth=1.0)\n";
    s << "s.unsetPrimaryObject()\n";

    s << "v1 = p2.vertices\n";
    s << "p2.ReferencePoint(point=v1[2])\n";


    // assembly
    s << "a1 = mdb.models['Model-1'].rootAssembly\n";
    s << "a1.DatumCsysByDefault(CARTESIAN)\n";

    // add and rotate main part
    s << "inst1 = a1.Instance(name='MyPart1-1', part=p, dependent=ON)\n";
    s << "a1.rotate(instanceList=('MyPart1-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=-90.0)\n";

    // add and rotate indenter
    s << "a1.Instance(name='Part-2-1', part=p2, dependent=ON)\n";
    s << "a1.translate(instanceList=('Part-2-1', ), vector=(0.0, 1.05, 0.0))\n";


    // create step
    s << "mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=2.0, improvedDtMethod=ON)\n";

    s << "mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=300)\n";

    // gravity load
    s << "mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Step-1',comp2=-10.0, distributionType=UNIFORM, field='')\n";

    // BC - pinned nodes
    s << "region = inst1.sets['Set3-pinned']\n";
    s << "mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', region=region, localCsys=None)\n";

    // BC - moving indenter
    s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
    s << "refPoints1=(r1[2], )\n";
    s << "region = a1.Set(referencePoints=refPoints1, name='Set-1-indenterRP')\n";
    s << "mdb.models['Model-1'].VelocityBC(name='BC-2', createStepName='Step-1', "
        "region=region, v1=0.0, v2=-0.001, v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, "
        "amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')\n";

    // rigid body constraint
    s << "s1 = a1.instances['Part-2-1'].faces\n";
    s << "side2Faces1 = s1[0:1]\n";
    s << "region5=a1.Surface(side2Faces=side2Faces1, name='Surf-1')\n";
    s << "r1 = a1.instances['Part-2-1'].referencePoints\n";
    s << "refPoints1=(r1[2], )\n";
    s << "region1=regionToolset.Region(referencePoints=refPoints1)\n";
    s << "mdb.models['Model-1'].RigidBody(name='Constraint-1', refPointRegion=region1, surfaceRegion=region5)\n";

    // create interaction property
    s << "mdb.models['Model-1'].ContactProperty('IntProp-1')\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior("
        "formulation=FRICTIONLESS)\n";
    s << "mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior("
        "pressureOverclosure=HARD, allowSeparation=ON, "
        "constraintEnforcementMethod=DEFAULT)\n";

    // create interaction itself
    s << "mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Step-1')\n";
    s << "mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep("
        "stepName='Step-1', useAllstar=ON)\n";
    s << "mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep("
        "stepName='Step-1', assignments=((GLOBAL, SELF, 'IntProp-1'), ))\n";


    // record indenter force

    s << "mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name="
        "'H-Output-2', rebar=EXCLUDE, region="
        "mdb.models['Model-1'].rootAssembly.sets['Set-1-indenterRP'], sectionPoints="
        "DEFAULT, timeInterval=0.0001, variables=('RF2', ))\n";

    //create job
    s << "mdb.Job(name='Job-1a', model='Model-1', description='', type=ANALYSIS,"
    "atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,"
    "memoryUnits=PERCENTAGE, explicitPrecision=SINGLE,"
    "nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF,"
    "contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',"
    "resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=4,"
    "activateLoadBalancing=False, numThreadsPerMpiProcess=1,"
    "multiprocessingMode=DEFAULT, numCpus=4)\n";

    s.close();
}


void icy::Mesh::RotateSample(double angleInDegrees)
{
    double alpha = angleInDegrees*M_PI/180.;
    double cosA = std::cos(alpha);
    double sinA = std::sin(alpha);
    Eigen::Matrix3d R;
    R << cosA, -sinA, 0,
            sinA, cosA, 0,
            0, 0, 1;
    for(Node *nd : nodes)
    {
        nd->x0 = nd->xn = nd->xt = R*nd->x0;
    }
}

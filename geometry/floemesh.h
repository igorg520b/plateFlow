#ifndef FL333_H
#define FL333_H

#include <vector>
#include <string>
#include <algorithm>

#include <spdlog/spdlog.h>

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>

#include <gmsh.h>
#include <Eigen/Core>
#include <H5Cpp.h>

#include <mutex>

#include "ConcurrentPool.h"
#include "equationofmotionsolver.h"



namespace icy { class FloeMesh; struct Node; struct Element; }

class icy::FloeMesh
{
public:
    FloeMesh();
    ~FloeMesh();
    FloeMesh& operator=(FloeMesh&) = delete;

    std::vector<icy::Node*> nodes;
    std::vector<icy::Element*> elems;
    std::mutex mtxMeshUpdate;

    void Reset();
    void GoToStep0();

    icy::Node* AddNode();
    icy::Element* AddElement();

private:
    constexpr static unsigned reserveConst = 100000;
    static ConcurrentPool<Node> NodeFactory;
    static ConcurrentPool<Element> ElementFactory;

};
#endif

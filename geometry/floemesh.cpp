#include "floemesh.h"
#include "element.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>



icy::ConcurrentPool<icy::Node> icy::FloeMesh::NodeFactory(reserveConst);
icy::ConcurrentPool<icy::Element> icy::FloeMesh::ElementFactory(reserveConst);


icy::FloeMesh::FloeMesh()
{
    nodes.reserve(reserveConst);
    elems.reserve(reserveConst);
}

icy::FloeMesh::~FloeMesh()
{
    Reset();
}

void icy::FloeMesh::Reset()
{
    NodeFactory.release(nodes);
    ElementFactory.release(elems);
}

void icy::FloeMesh::GoToStep0()
{
    for(Node *nd : nodes)
    {
        nd->xn = nd->x0;
        nd->vn.setZero();
    }
}


icy::Node* icy::FloeMesh::AddNode()
{
    Node* nd = NodeFactory.take();
    nd->Reset();
    nd->globId = (int)nodes.size();
    nodes.push_back(nd);
    return nd;
}

icy::Element* icy::FloeMesh::AddElement()
{
    Element* elem = ElementFactory.take();
    elem->Reset();
    elems.push_back(elem);
    return elem;
}


void icy::FloeMesh::Import(std::string fileName)
{

}

void icy::FloeMesh::HDF5LoadStep(H5::H5File *file, int step)
{

}

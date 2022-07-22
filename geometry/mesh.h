#ifndef MESH_H
#define MESH_H


#include <gmsh.h>

#include <vector>

#include <Eigen/Core>


namespace icy {class Mesh;}

class icy::Mesh
{
public:
    Mesh();
    ~Mesh();
    Mesh& operator=(Mesh&) = delete;

    void Reset();

};

#endif // MESH_H

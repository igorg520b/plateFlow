#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing tbb headers
#ifndef FL333_H
#define FL333_H

#include <vector>
#include <string>
#include <algorithm>

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>

#include <gmsh.h>
#include <Eigen/Core>
#include <H5Cpp.h>

#include <QMutex>

#include "ConcurrentPool.h"
#include "equationofmotionsolver.h"

namespace icy { class Mesh; struct Node; struct Element; }

class icy::Mesh
{
public:
    Mesh();
    ~Mesh();
    Mesh& operator=(Mesh&) = delete;


};
#endif
#endif //Q_MOC_RUN

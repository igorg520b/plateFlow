#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
#ifndef NODE_H
#define NODE_H

#include <functional>
#include <vector>
#include <Eigen/Core>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"

namespace icy { struct Node; struct Element; struct CohesiveZone; class MeshFragment;}

struct icy::Node
{
    Node() { Reset(); };
    ~Node() = default;
    Node& operator=(Node&) = delete;

    void Reset();
    void Initialize(double x, double y);
    void Initialize(const Node *other);

};

#endif // NODE_H
#endif // Q_MOC_RUN

#ifndef NODE_H
#define NODE_H

#include <functional>
#include <vector>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"
#include <Eigen/Core>

namespace icy { struct Node; struct Element; struct CohesiveZone; class MeshFragment;}

struct icy::Node
{
    Node() { Reset(); };
    ~Node() = default;
    Node& operator=(Node&) = delete;

    void Reset();
    void Initialize(double x, double y);
    void Initialize(const Node *other);

    int eqId;      // squential number in the system of equations (-1 if prescribed)
    int globId;
    int locId;          // sequential number in a given floe
    double area;        // mass that the node "represents", for applying various forces
    double vertical_force; // for testing
    bool isBoundary;

    // initial configuration
    Eigen::Matrix<double,5,1> x0;

    // at step n: displacement, position, velocity, acceleration
    Eigen::Matrix<double,5,1> un, xn, vn, an;
    Eigen::Matrix<double,5,1> ut, xt, vt, at; // at step n+1

};

#endif // NODE_H

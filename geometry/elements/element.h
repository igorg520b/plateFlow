#ifndef ELEMENT123_H
#define ELEMENT123_H

#include "node.h"
#include "equationofmotionsolver.h"
#include "parameters_sim.h"

namespace icy { struct Element; struct Node; }

struct icy::Element
{
    void Reset();
    void Precompute();
    void AddToSparsityStructure(EquationOfMotionSolver &eq) const;
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, const Params &prms, const double timeStep);
    void ComputeVisualizedVariables(const Params &prms);

};

#endif // ELEMENT123_H

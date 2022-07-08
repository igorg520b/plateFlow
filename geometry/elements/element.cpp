#include "element.h"
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <spdlog/spdlog.h>
#include <complex>
#include <cmath>

void icy::Element::Reset()
{}

void icy::Element::Precompute()
{}

void icy::Element::AddToSparsityStructure(EquationOfMotionSolver &eq) const
{}

bool icy::Element::ComputeEquationEntries(EquationOfMotionSolver &eq, const Params &prms, const double timeStep)
{}

void icy::Element::ComputeVisualizedVariables(const Params &prms)
{}


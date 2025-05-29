#pragma once
#include "NonlinearSystem.h"
#include "NumericalJacobian.h"
#include "LinearAlgebra.h"

class LevenbergMarquardtSolver {
    const NonlinearSystem& system_;
    NumericalJacobian jacobian_;
public:
    LevenbergMarquardtSolver(const NonlinearSystem& sys);
    bool solve(Vec& x);
};
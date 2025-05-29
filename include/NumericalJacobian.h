#pragma once
#include "NonlinearSystem.h"

constexpr double Tol = 1e-8;
using Mat = double[N][N];

struct NumericalJacobian {
    static constexpr double h = Tol;
    void compute(const NonlinearSystem& sys, const Vec& x, Mat& jacobian) const;
};
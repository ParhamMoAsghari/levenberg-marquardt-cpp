#pragma once
constexpr int N = 5;
using Vec = double[N];

struct NonlinearSystem {
    void computeF(const Vec& x, Vec& residual) const;
};
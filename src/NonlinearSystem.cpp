#include "NonlinearSystem.h"
#include <cmath>

void NonlinearSystem::computeF(const Vec& x, Vec& residual) const {
    residual[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] - 10;
    residual[1] = x[0]*x[1] + x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4]*x[0] - 5;
    residual[2] = std::exp(x[0]) + std::log(x[1] + 1.1) - x[2]*x[2] + std::sin(x[3]) - x[4];
    residual[3] = x[0]*x[4] - std::cos(x[1]) + x[2]*x[3] - 1;
    residual[4] = x[0] + x[1]*x[1] - x[2] + x[3]*x[3]*x[3] - std::exp(-x[4]) - 2;
}
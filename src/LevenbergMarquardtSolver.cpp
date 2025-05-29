#include "LevenbergMarquardtSolver.h"
#include <iostream>
#include <cstring>

constexpr int MAX_ITER = 1000;
constexpr double TOL = 1e-6;

LevenbergMarquardtSolver::LevenbergMarquardtSolver(const NonlinearSystem& sys) : system_(sys) {}

bool LevenbergMarquardtSolver::solve(Vec& x) {
    double lambda = 0.01;
    const double v = 2.0;
    Vec residual, residual_new, step, x_new, gradient;
    Mat J, Jt, A, A_lm;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        system_.computeF(x, residual);
        double err = norm(residual);
        if (err < TOL) {
            std::cout << "Converged at iteration " << iter << ", error: " << err << "\n";
            return true;
        }

        jacobian_.compute(system_, x, J);
        transpose(J, Jt);

        // compute gradient: g = Jt * f
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j)
                sum += Jt[i][j] * residual[j];
            gradient[i] = -sum;
        }

        matMul(Jt, J, A);
        std::memcpy(A_lm, A, sizeof(A));
        addLambdaIdentity(A_lm, lambda);

        if (!gaussSolve(A_lm, gradient, step)) {
            std::cout << "Failed linear solve\n";
            return false;
        }

        for (int i = 0; i < N; ++i)
            x_new[i] = x[i] + step[i];

        system_.computeF(x_new, residual_new);
        double err_new = norm(residual_new);

        if (err_new < err) {
            std::memcpy(x, x_new, sizeof(Vec));
            lambda /= v;
            if (lambda < 1e-7) lambda = 1e-7;
        } else {
            lambda *= v;
            if (lambda > 1e7) {
                std::cout << "Lambda too large, stopping\n";
                return false;
            }
        }
    }
    std::cout << "Max iterations reached, final error: " << norm(residual) << "\n";
    return false;
}
#include <iostream>
#include "NonlinearSystem.h"
#include "LevenbergMarquardtSolver.h"

int main() {
    NonlinearSystem sys;
    LevenbergMarquardtSolver solver(sys);

    Vec x = {1, 1, 1, 1, 1};

    if (solver.solve(x)) {
        std::cout << "Solution:\n";
        for (int i = 0; i < N; ++i)
            std::cout << "x[" << i << "] = " << x[i] << "\n";
    } else {
        std::cout << "No convergence\n";
    }
    return 0;
}
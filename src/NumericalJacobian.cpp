#include "NumericalJacobian.h"
#include <cstring>

void NumericalJacobian::compute(const NonlinearSystem& sys, const Vec& x, Mat& jacobian) const {
    Vec f0, f1;
    sys.computeF(x, f0);

    for (int j = 0; j < N; ++j) {
        Vec xh;
        std::memcpy(xh, x, sizeof(Vec));
        xh[j] += h;
        sys.computeF(xh, f1);
        for (int i = 0; i < N; ++i) {
            jacobian[i][j] = (f1[i] - f0[i]) / h;
        }
    }
}
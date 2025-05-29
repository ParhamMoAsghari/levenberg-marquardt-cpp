#include "LinearAlgebra.h"

void matVecMul(const Mat& A, const Vec& x, Vec& y) {
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        const double* Ai = A[i];
        for (int j = 0; j < N; ++j) {
            sum += Ai[j] * x[j];
        }
        y[i] = sum;
    }
}

void matMul(const Mat& A, const Mat& B, Mat& C) {
    for (int i = 0; i < N; ++i) {
        double* Ci = C[i];
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k)
                sum += A[i][k] * B[k][j];
            Ci[j] = sum;
        }
    }
}

void transpose(const Mat& A, Mat& At) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            At[j][i] = A[i][j];
}

void addLambdaIdentity(Mat& A, double lambda) {
    for (int i = 0; i < N; ++i)
        A[i][i] += lambda;
}

bool gaussSolve(Mat A, Vec b, Vec& x) {
    for (int i = 0; i < N; ++i) {
        int pivot = i;
        double max_val = std::abs(A[i][i]);
        for (int k = i + 1; k < N; ++k) {
            double val = std::abs(A[k][i]);
            if (val > max_val) {
                max_val = val;
                pivot = k;
            }
        }
        if (max_val < 1e-14) return false;
        if (pivot != i) {
            for (int j = 0; j < N; ++j) std::swap(A[i][j], A[pivot][j]);
            std::swap(b[i], b[pivot]);
        }

        double Ai_i = A[i][i];
        for (int k = i + 1; k < N; ++k) {
            double f = A[k][i] / Ai_i;
            for (int j = i; j < N; ++j) {
                A[k][j] -= f * A[i][j];
            }
            b[k] -= f * b[i];
        }
    }

    for (int i = N - 1; i >= 0; --i) {
        double Ai_i = A[i][i];
        if (std::abs(Ai_i) < 1e-14) return false;
        double sum = b[i];
        for (int j = i + 1; j < N; ++j)
            sum -= A[i][j] * x[j];
        x[i] = sum / Ai_i;
    }
    return true;
}

double norm(const Vec& v) {
    double sum = 0.0;
    for (int i = 0; i < N; ++i)
        sum += v[i] * v[i];
    return std::sqrt(sum);
}
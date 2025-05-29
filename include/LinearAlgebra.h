#pragma once
#include <cmath>
#include <cstring>
#include <utility>
#include "NumericalJacobian.h"

void matVecMul(const Mat& A, const Vec& x, Vec& y);
void matMul(const Mat& A, const Mat& B, Mat& C);
void transpose(const Mat& A, Mat& At);
void addLambdaIdentity(Mat& A, double lambda);
bool gaussSolve(Mat A, Vec b, Vec& x);
double norm(const Vec& v);
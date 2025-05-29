# Levenberg-Marquardt Algorithm Implementation

This project implements the Levenberg-Marquardt nonlinear least squares solver in C++. It is designed with modular code for clarity and extensibility.

## Features

- Matrix and vector algebra utilities
- Numerical Jacobian calculation
- Levenberg-Marquardt optimization solver
- Example usage in `main.cpp`

## Build Instructions

Requires CMake and a C++17 compatible compiler.

```bash
mkdir build
cd build
cmake ..
cmake --build .
./lm_solver

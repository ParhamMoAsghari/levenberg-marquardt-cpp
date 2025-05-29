cmake_minimum_required(VERSION 3.10)
project(LevenbergMarquardtSolver LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# Include headers
include_directories(${CMAKE_SOURCE_DIR}/include)

# Collect sources
set(SRC
    src/NonlinearSystem.cpp
    src/NumericalJacobian.cpp
    src/LinearAlgebra.cpp
    src/LevenbergMarquardtSolver.cpp
    main.cpp
)

# Define executable
add_executable(lm_solver ${SRC})

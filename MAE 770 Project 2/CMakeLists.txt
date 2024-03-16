cmake_minimum_required(VERSION 3.10)
project("MAE 770 Project 2")

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags for Release builds
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")

# Eigen3 include directory
set(EIGEN3_INCLUDE_DIR C:/Users/luis2/eigen-3.4.0/eigen-3.4.0)

# Include directories for your project and Eigen3
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Create an executable with all source files
add_executable(NOZZLE src/main.cpp  "src/WriteData.cpp" "src/MeshProccessing.cpp" "src/SpeciesThermo.cpp" "src/CellStateVariables.cpp" "src/SimParameters.cpp" "src/BuildCellJacobians.cpp")

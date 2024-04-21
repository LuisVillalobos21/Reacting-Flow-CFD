cmake_minimum_required(VERSION 3.10)
project("MAE 770 Project 2")

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Append /O2 and /DEBUG for Release builds without overwriting existing flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox" CACHE STRING "Release flags" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG" CACHE STRING "Linker flags for Release" FORCE)

# Eigen3 include directory
set(EIGEN3_INCLUDE_DIR "C:/Users/luis2/eigen-3.4.0/eigen-3.4.0")

# Include directories for your project and Eigen3
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Create an executable with all source files
add_executable(SOLVER 
  src/main.cpp  
  src/WriteData.cpp 
  src/MeshProccessing.cpp 
  src/SpeciesThermo.cpp 
  src/CellStateVariables.cpp 
  src/SimParameters.cpp 
  src/BuildCellJacobians.cpp 
  src/RightHandSide.cpp 
  src/TimeIntegration.cpp
 "src/Chemistry.cpp")

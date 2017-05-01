project(UnscentedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++0x)

set(sources
   src/ukf.cpp
   src/main.cpp
   src/tools.cpp
   src/cmaes/cmaes.h
   src/cmaes/cmaes_interface.h
   src/cmaes/cmaes.cpp)

add_executable(UnscentedKF ${sources})

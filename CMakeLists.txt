cmake_minimum_required(VERSION 3.25)
project(examen_JM)

set(CMAKE_CXX_STANDARD 17)

include(cmake/CPM.cmake)

#https://github.com/fmtlib/fmt.git
CPMAddPackage("gh:fmtlib/fmt#10.0.0")

#https://github.com/abdes/cryptopp-cmake.git
CPMAddPackage("gh:abdes/cryptopp-cmake#v8.8.0")
find_package(MPI REQUIRED)

add_executable(examen_JM main.cpp)
target_link_libraries(examen_JM fmt::fmt MPI::MPI_CXX)

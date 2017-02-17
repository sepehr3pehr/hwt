cmake_minimum_required(VERSION 2.6)
project (nbt)

find_library(HDF5_LIB hdf5)

set (CMAKE_CXX_FLAGS "-DSINGLE_CORE -O2 -fopenmp -DNDEBUG -std=c++0x")

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories("/usr/include/hdf5/serial/")
add_executable (nbt src/node.cpp src/root_node.cpp src/loadVar.cpp test/test_create_tree.cpp)
target_link_libraries (nbt ${HDF5_LIB})



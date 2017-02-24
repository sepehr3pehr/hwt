cmake_minimum_required(VERSION 2.6)
project (hwt)

find_library(HDF5_LIB hdf5)

set (CMAKE_CXX_FLAGS "-DSINGLE_CORE -O2 -fopenmp -DNDEBUG -std=c++0x")

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories("/usr/include/hdf5/serial/")


add_executable (hwt src/node.cpp src/root_node.cpp src/loadVar.cpp interface/hwt_interface.cpp)
target_link_libraries (hwt ${HDF5_LIB})

project (linscan)
add_executable (linscan interface/linscan_interface.cpp src/linscan.cpp src/loadVar.cpp)
target_link_libraries (linscan ${HDF5_LIB})


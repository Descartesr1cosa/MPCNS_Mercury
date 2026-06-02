cmake_minimum_required(VERSION 3.5)
project(Z1_Topo_Test)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories(../../include)
include_directories(.)

aux_source_directory(../../src/0_basic DIRSRCS0)
aux_source_directory(../../src/1_grid DIRSRCS1)
aux_source_directory(../../src/2_topology DIRSRCS2)

add_executable(topo_test
    topo_test.cpp
    ${DIRSRCS0}
    ${DIRSRCS1}
    ${DIRSRCS2})

target_link_libraries(topo_test PRIVATE ${MPI_LIBRARIES} stdc++fs)

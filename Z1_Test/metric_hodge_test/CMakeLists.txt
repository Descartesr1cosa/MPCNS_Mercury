cmake_minimum_required(VERSION 3.5)
project(Z1_Metric_Hodge_Test)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories(../../include)
include_directories(../common)
include_directories(.)

aux_source_directory(../../src/0_basic DIRSRCS0)
aux_source_directory(../../src/1_grid DIRSRCS1)
aux_source_directory(../../src/2_topology DIRSRCS2)
aux_source_directory(../../src/3_field DIRSRCS3)
aux_source_directory(../../src/7_metric DIRSRCS7)

add_executable(metric_hodge_test
    metric_hodge_test.cpp
    ${DIRSRCS0}
    ${DIRSRCS1}
    ${DIRSRCS2}
    ${DIRSRCS3}
    ${DIRSRCS7})

target_link_libraries(metric_hodge_test PRIVATE ${MPI_LIBRARIES} stdc++fs)

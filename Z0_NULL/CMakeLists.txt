cmake_minimum_required(VERSION 3.5)
project(Z0_NULL)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json" FORCE)

find_package(MPI REQUIRED)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

include_directories(${MPI_INCLUDE_PATH})
include_directories(../include)
include_directories(./include)
include_directories(.)

aux_source_directory(../src/0_basic DIRSRCS0)
aux_source_directory(../src/1_grid DIRSRCS1)
aux_source_directory(../src/2_topology DIRSRCS2)
aux_source_directory(../src/3_field DIRSRCS3)
aux_source_directory(../src/4_halo DIRSRCS4)
aux_source_directory(../src/7_metric DIRSRCS7)
aux_source_directory(./src Z0_NULL_SRCS)
aux_source_directory(./src/boundary Z0_NULL_BOUNDARY_SRCS)
aux_source_directory(./src/physics Z0_NULL_PHYSICS_SRCS)
add_executable(z0_null
    main.cpp
    ${DIRSRCS0}
    ${DIRSRCS1}
    ${DIRSRCS2}
    ${DIRSRCS3}
    ${DIRSRCS4}
    ${DIRSRCS7}
    ${Z0_NULL_SRCS}
    ${Z0_NULL_BOUNDARY_SRCS}
    ${Z0_NULL_PHYSICS_SRCS})

target_link_libraries(z0_null PRIVATE ${MPI_LIBRARIES} stdc++fs)

cmake_minimum_required(VERSION 3.5)
project(Z1_NULL)

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
add_executable(z1_null
    main.cpp
    ./src/Z1_Control.cpp
    ${DIRSRCS0}
    ${DIRSRCS1}
    ${DIRSRCS2})

target_link_libraries(z1_null PRIVATE ${MPI_LIBRARIES} stdc++fs)

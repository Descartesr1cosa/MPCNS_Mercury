cmake_minimum_required(VERSION 3.5)
project(MPCNS_TEST)
set(CMAKE_CXX_STANDARD 20)

# Find MPI
find_package(MPI REQUIRED)

# Set Compile Options
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# include files directory
include_directories(${MPI_INCLUDE_PATH})
include_directories(../include)
include_directories(.)

# source files directiory
aux_source_directory(../src/0_basic  DIRSRCS0)
aux_source_directory(../src/1_grid  DIRSRCS1)
aux_source_directory(../src/2_topology  DIRSRCS2)
aux_source_directory(../src/3_field  DIRSRCS3)

# build the executable program
add_executable(MPCNS_TEST main.cpp ${DIRSRCS0} ${DIRSRCS1} ${DIRSRCS2} ${DIRSRCS3})

# set compile options
#target_compile_options(MPCNS_TEST PRIVATE -std=c++20 -fconcepts)

# add the mpi linkings
target_link_libraries(MPCNS_TEST ${MPI_LIBRARIES} stdc++fs)

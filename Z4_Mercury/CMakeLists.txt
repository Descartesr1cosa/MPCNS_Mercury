cmake_minimum_required(VERSION 3.5)
project(MercuryZ4)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json" FORCE)

# Find MPI
find_package(MPI REQUIRED)

# Find PETSc
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

# Set Compile Options
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# include files directory
include_directories(${MPI_INCLUDE_PATH})
include_directories(../include)
include_directories(../Z4_Mercury/include)
include_directories(.)

# source files directiory
# Basic Files
aux_source_directory(../src/0_basic  DIRSRCS0)
aux_source_directory(../src/1_grid  DIRSRCS1)
aux_source_directory(../src/2_topology  DIRSRCS2)
aux_source_directory(../src/3_field  DIRSRCS3)
aux_source_directory(../src/4_halo  DIRSRCS4)
aux_source_directory(../src/5_boundary  DIRSRCS5)
aux_source_directory(../src/6_io  DIRSRCS6)
aux_source_directory(../src/7_metric  DIRSRCS7)

# Mercury Solver Files
aux_source_directory(./src  SRC)
aux_source_directory(./src/boundary  SRC_BOUND)
aux_source_directory(./src/operators  SRC_OPERATOR)
aux_source_directory(./src/physics  SRC_PHY)
# aux_source_directory(./src/3_io  SRC2)
# aux_source_directory(./src/4_solver  SRC3)
aux_source_directory(./src/tools  SRC_TOOL)
aux_source_directory(./src/implicit  SRC_IMP)
set(SRC_SETUP ./src/MercurySolver_Setup.cpp)

# build the executable program
add_executable(MercuryZ4 main.cpp 
${DIRSRCS0} ${DIRSRCS1} ${DIRSRCS2} ${DIRSRCS3}
 ${DIRSRCS4} ${DIRSRCS5} ${DIRSRCS6} ${DIRSRCS7} ${SRC}
 ${SRC_OPERATOR} ${SRC_BOUND} ${SRC_PHY} ${SRC_TOOL} ${SRC_IMP} ${SRC_SETUP})
#  ${SRC0} ${SRC1} ${SRC2} ${SRC3} ${SRC4} ${SRC5} ${SRC6})

# add the mpi linkings
target_link_libraries(MercuryZ4 PRIVATE ${MPI_LIBRARIES} PkgConfig::PETSC stdc++fs)

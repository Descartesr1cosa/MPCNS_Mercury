cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0069 NEW)
project(MercuryZ4)
option(MERCURY_SANITIZE_ADDRESS "Build MercuryZ4 with AddressSanitizer" OFF)
option(OUTPUT_DEC_JEDGE "Allow debug-only solver-native DEC Jedge restart fields" OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json" FORCE)

# Find MPI
find_package(MPI REQUIRED COMPONENTS CXX)

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

# Shared infrastructure sources.
aux_source_directory(../src/0_basic  DIRSRCS0)
aux_source_directory(../src/1_grid  DIRSRCS1)
aux_source_directory(../src/2_topology  DIRSRCS2)
aux_source_directory(../src/3_field  DIRSRCS3)
aux_source_directory(../src/4_halo/core DIRSRCS4_CORE)
aux_source_directory(../src/4_halo/build DIRSRCS4_BUILD)
aux_source_directory(../src/4_halo/component_copy DIRSRCS4_COMPONENT_COPY)
aux_source_directory(../src/4_halo/forms DIRSRCS4_FORMS)
aux_source_directory(../src/4_halo/coupling DIRSRCS4_COUPLING)
aux_source_directory(../src/4_halo/owner DIRSRCS4_OWNER)
aux_source_directory(../src/5_boundary  DIRSRCS5)
aux_source_directory(../src/6_io  DIRSRCS6)
aux_source_directory(../src/7_metric  DIRSRCS7)
aux_source_directory(../src/8_dec  DIRSRCS8)

# Mercury solver sources.
aux_source_directory(./src  SRC)
aux_source_directory(./src/boundary  SRC_BOUND)
aux_source_directory(./src/operators  SRC_OPERATOR)
aux_source_directory(./src/physics  SRC_PHY)
aux_source_directory(./src/tools  SRC_TOOL)
set(SRC_IMP ./src/implicit/Mercury_ImplicitResistive.cpp)

# Build the executable.
add_executable(MercuryZ4 main.cpp 
${DIRSRCS0} ${DIRSRCS1} ${DIRSRCS2} ${DIRSRCS3}
 ${DIRSRCS4_CORE} ${DIRSRCS4_BUILD} ${DIRSRCS4_COMPONENT_COPY}
 ${DIRSRCS4_FORMS} ${DIRSRCS4_COUPLING} ${DIRSRCS4_OWNER}
 ${DIRSRCS5} ${DIRSRCS6} ${DIRSRCS7} ${DIRSRCS8} ${SRC}
 ${SRC_OPERATOR} ${SRC_BOUND} ${SRC_PHY} ${SRC_TOOL} ${SRC_IMP})

# Link MPI and PETSc.
target_link_libraries(MercuryZ4 PRIVATE ${MPI_LIBRARIES} PkgConfig::PETSC stdc++fs)
if(OUTPUT_DEC_JEDGE)
    target_compile_definitions(MercuryZ4 PRIVATE OUTPUT_DEC_JEDGE=1)
endif()

if(MERCURY_SANITIZE_ADDRESS)
    target_compile_options(MercuryZ4 PRIVATE -O1 -g -fno-omit-frame-pointer -fsanitize=address)
    target_link_options(MercuryZ4 PRIVATE -fsanitize=address)
endif()

# Allow cross-translation-unit inlining of the small Field/FieldBlock accessors
# used throughout the inner loops.  This changes code generation only; the
# numerical algorithm and floating-point model remain unchanged.
include(CheckIPOSupported)
check_ipo_supported(RESULT MERCURY_IPO_SUPPORTED OUTPUT MERCURY_IPO_ERROR)
if(MERCURY_IPO_SUPPORTED AND NOT MERCURY_SANITIZE_ADDRESS)
    set_property(TARGET MercuryZ4 PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

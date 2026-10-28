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
# Keep the initialization-only prepared metadata builder explicit: source
# discovery alone does not make a newly added file trigger CMake regeneration.
list(APPEND DIRSRCS4_CORE ../src/4_halo/core/HaloPreparedMetadata.cpp)
list(APPEND DIRSRCS4_CORE ../src/4_halo/core/HaloPreparedExchange.cpp)
list(REMOVE_DUPLICATES DIRSRCS4_CORE)
aux_source_directory(../src/4_halo/build DIRSRCS4_BUILD)
aux_source_directory(../src/4_halo/component_copy DIRSRCS4_COMPONENT_COPY)
aux_source_directory(../src/4_halo/forms DIRSRCS4_FORMS)
aux_source_directory(../src/4_halo/coupling DIRSRCS4_COUPLING)
list(APPEND DIRSRCS4_COUPLING ../src/4_halo/coupling/HaloCouplingPrepared.cpp)
aux_source_directory(../src/4_halo/owner DIRSRCS4_OWNER)
aux_source_directory(../src/5_boundary  DIRSRCS5)
aux_source_directory(../src/6_io  DIRSRCS6)
aux_source_directory(../src/7_metric  DIRSRCS7)
aux_source_directory(../src/8_dec  DIRSRCS8)
set(DIRSRCS9 ../src/9_runtime/SimulationContext.cpp)

if(NOT TARGET mpcns_mesh)
    add_library(mpcns_mesh STATIC ../src/mpcns/mesh/mesh_partition.cpp)
    target_include_directories(mpcns_mesh PUBLIC ../include)
endif()

if(NOT TARGET mpcns_parallel)
    add_library(mpcns_parallel STATIC
        ../src/mpcns/parallel/exchange_plan.cpp
        ../src/mpcns/parallel/oriented_form_plan.cpp
        ../src/mpcns/parallel/oriented_corner_plan.cpp
        ../src/mpcns/parallel/oriented_vertex_plan.cpp
        ../src/mpcns/parallel/face_2form_edge_plan.cpp
        ../src/mpcns/parallel/face_2form_vertex_plan.cpp)
    target_include_directories(mpcns_parallel PUBLIC ../include)
    target_link_libraries(mpcns_parallel PUBLIC mpcns_mesh)
endif()

if(NOT TARGET mpcns_legacy_mesh)
    add_library(mpcns_legacy_mesh STATIC
        ../src/mpcns/compat/legacy_mesh_partition_adapter.cpp
        ../src/mpcns/compat/legacy_halo_plan_adapter.cpp
        ../src/mpcns/compat/legacy_oriented_form_plan_adapter.cpp
        ../src/mpcns/compat/legacy_oriented_corner_plan_adapter.cpp
        ../src/mpcns/compat/legacy_oriented_vertex_plan_adapter.cpp
        ../src/mpcns/compat/legacy_face_2form_edge_plan_adapter.cpp
        ../src/mpcns/compat/legacy_face_2form_vertex_plan_adapter.cpp)
    target_include_directories(mpcns_legacy_mesh PUBLIC ../include)
    target_link_libraries(mpcns_legacy_mesh PUBLIC mpcns_mesh mpcns_parallel)
endif()

if(NOT TARGET mpcns_state)
    add_library(mpcns_state STATIC
        ../src/mpcns/state/state_layout.cpp
        ../src/mpcns/state/simulation_state.cpp)
    target_include_directories(mpcns_state PUBLIC ../include)
endif()

if(NOT TARGET mpcns_state_compat)
    add_library(mpcns_state_compat STATIC
        ../src/mpcns/compat/legacy_state_layout.cpp
        ../src/mpcns/compat/legacy_state_view_adapter.cpp)
    target_include_directories(mpcns_state_compat PUBLIC ../include)
    target_link_libraries(mpcns_state_compat PUBLIC mpcns_state)
endif()

# New bottom-layer modules are explicit targets rather than source globs.
if(NOT TARGET mpcns_discrete)
    add_library(mpcns_discrete STATIC
        ../src/mpcns/discrete/incidence_operator.cpp
        ../src/mpcns/discrete/lumped_hodge.cpp
        ../src/mpcns/discrete/dual_face_inverse_hodge.cpp)
    target_include_directories(mpcns_discrete PUBLIC ../include)
    target_link_libraries(mpcns_discrete PUBLIC mpcns_state)
endif()

if(NOT TARGET mpcns_operators)
    add_library(mpcns_operators STATIC
        ../src/mpcns/operators/electromagnetic_operators.cpp
        ../src/mpcns/operators/faraday_operator.cpp
        ../src/mpcns/operators/reconstruction_operators.cpp
        ../src/mpcns/operators/pole_current_constraint.cpp
        ../src/mpcns/constraints/least_squares_3.cpp)
    target_include_directories(mpcns_operators PUBLIC ../include)
    target_link_libraries(mpcns_operators PUBLIC mpcns_discrete)
endif()

if(NOT TARGET mpcns_geometry)
    add_library(mpcns_geometry STATIC
        ../src/mpcns/geometry/geometry_primitives.cpp
        ../src/mpcns/geometry/geometry_snapshot.cpp
        ../src/mpcns/geometry/geometry_builder.cpp
        ../src/mpcns/geometry/singular_edge_geometry.cpp
        ../src/mpcns/geometry/reconstruction_kernels.cpp)
    target_include_directories(mpcns_geometry PUBLIC ../include)
endif()

if(NOT TARGET mpcns_legacy_geometry)
    add_library(mpcns_legacy_geometry STATIC
        ../src/mpcns/compat/legacy_geometry_field_ids.cpp
        ../src/mpcns/compat/legacy_geometry_field_registration.cpp
        ../src/mpcns/compat/legacy_geometry_snapshot.cpp
        ../src/mpcns/compat/legacy_geometry_build_pipeline.cpp
        ../src/mpcns/compat/legacy_lumped_metric_publisher.cpp
        ../src/mpcns/compat/legacy_primal_metric_publisher.cpp
        ../src/mpcns/compat/legacy_reconstruction_metric_publisher.cpp
        ../src/mpcns/compat/legacy_field_reconstruction.cpp
        ../src/mpcns/compat/legacy_metric_diagnostics.cpp
        ../src/mpcns/compat/legacy_singular_edge_metric_publisher.cpp
        ../src/mpcns/compat/legacy_singular_edge_collector.cpp)

    target_include_directories(mpcns_legacy_geometry PUBLIC ../include)
    target_link_libraries(mpcns_legacy_geometry PUBLIC mpcns_discrete mpcns_geometry)
endif()

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
 ${DIRSRCS5} ${DIRSRCS6} ${DIRSRCS7} ${DIRSRCS8} ${DIRSRCS9} ${SRC}
 ${SRC_OPERATOR} ${SRC_BOUND} ${SRC_PHY} ${SRC_TOOL} ${SRC_IMP})

# Link MPI and PETSc.
target_link_libraries(MercuryZ4 PRIVATE
    mpcns_legacy_mesh
    mpcns_discrete
    mpcns_state_compat
    mpcns_operators
    mpcns_legacy_geometry
    ${MPI_LIBRARIES}
    PkgConfig::PETSC
    stdc++fs)
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

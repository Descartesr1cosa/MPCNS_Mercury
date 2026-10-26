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
include_directories(./test/common)
include_directories(.)

# Z0 tests run from the build directory.  Keep their compact input case next
# to the executable so validation is reproducible without a manual copy step.
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/CASE" DESTINATION "${CMAKE_BINARY_DIR}")

aux_source_directory(../src/0_basic DIRSRCS0)
aux_source_directory(../src/1_grid DIRSRCS1)
aux_source_directory(../src/2_topology DIRSRCS2)
aux_source_directory(../src/3_field DIRSRCS3)
aux_source_directory(../src/4_halo/core DIRSRCS4_CORE)
aux_source_directory(../src/4_halo/build DIRSRCS4_BUILD)
aux_source_directory(../src/4_halo/component_copy DIRSRCS4_COMPONENT_COPY)
aux_source_directory(../src/4_halo/forms DIRSRCS4_FORMS)
aux_source_directory(../src/4_halo/coupling DIRSRCS4_COUPLING)
aux_source_directory(../src/4_halo/owner DIRSRCS4_OWNER)
aux_source_directory(../src/6_io DIRSRCS6)
aux_source_directory(../src/7_metric DIRSRCS7)
aux_source_directory(../src/8_dec DIRSRCS8)

if(NOT TARGET mpcns_discrete)
    add_library(mpcns_discrete STATIC
        ../src/mpcns/discrete/incidence_operator.cpp
        ../src/mpcns/discrete/lumped_hodge.cpp
        ../src/mpcns/discrete/dual_face_inverse_hodge.cpp)
    target_include_directories(mpcns_discrete PUBLIC ../include)
endif()

if(NOT TARGET mpcns_operators)
    add_library(mpcns_operators STATIC
        ../src/mpcns/operators/electromagnetic_operators.cpp
        ../src/mpcns/operators/faraday_operator.cpp)
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
aux_source_directory(./src Z0_NULL_SRCS)
aux_source_directory(./src/boundary Z0_NULL_BOUNDARY_SRCS)
aux_source_directory(./src/physics Z0_NULL_PHYSICS_SRCS)
aux_source_directory(./test/halo_comm Z0_NULL_TEST_HALO_SRCS)
aux_source_directory(./test/dec_chain Z0_NULL_TEST_DEC_SRCS)
aux_source_directory(./test/physical_io Z0_NULL_TEST_IO_SRCS)
aux_source_directory(./test/topology Z0_NULL_TEST_TOPOLOGY_SRCS)
add_executable(z0_null
    main.cpp
    ${DIRSRCS0}
    ${DIRSRCS1}
    ${DIRSRCS2}
    ${DIRSRCS3}
    ${DIRSRCS4_CORE}
    ${DIRSRCS4_BUILD}
    ${DIRSRCS4_COMPONENT_COPY}
    ${DIRSRCS4_FORMS}
    ${DIRSRCS4_COUPLING}
    ${DIRSRCS4_OWNER}
    ${DIRSRCS6}
    ${DIRSRCS7}
    ${DIRSRCS8}
    ${Z0_NULL_SRCS}
    ${Z0_NULL_BOUNDARY_SRCS}
    ${Z0_NULL_PHYSICS_SRCS}
    ${Z0_NULL_TEST_HALO_SRCS}
    ${Z0_NULL_TEST_DEC_SRCS}
    ${Z0_NULL_TEST_IO_SRCS}
    ${Z0_NULL_TEST_TOPOLOGY_SRCS})

target_link_libraries(z0_null PRIVATE
    mpcns_discrete
    mpcns_operators
    mpcns_legacy_geometry
    ${MPI_LIBRARIES}
    stdc++fs)

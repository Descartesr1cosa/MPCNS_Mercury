cmake_minimum_required(VERSION 3.16)
project(MarsZ5 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(MPI REQUIRED COMPONENTS CXX)

include_directories(${MPI_CXX_INCLUDE_DIRS})
include_directories(../include)
include_directories(include)

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
aux_source_directory(../src/5_boundary DIRSRCS5)
aux_source_directory(../src/6_io DIRSRCS6)
aux_source_directory(../src/7_metric DIRSRCS7)
aux_source_directory(../src/8_dec DIRSRCS8)
aux_source_directory(src MARS_SRC)

add_executable(MarsZ5 main.cpp
  ${DIRSRCS0} ${DIRSRCS1} ${DIRSRCS2} ${DIRSRCS3}
  ${DIRSRCS4_CORE} ${DIRSRCS4_BUILD} ${DIRSRCS4_COMPONENT_COPY}
  ${DIRSRCS4_FORMS} ${DIRSRCS4_COUPLING} ${DIRSRCS4_OWNER}
  ${DIRSRCS5} ${DIRSRCS6} ${DIRSRCS7} ${DIRSRCS8}
  ${MARS_SRC})

target_link_libraries(MarsZ5 PRIVATE MPI::MPI_CXX stdc++fs)

include(CheckIPOSupported)
check_ipo_supported(RESULT MARS_IPO_SUPPORTED OUTPUT MARS_IPO_ERROR)
if(MARS_IPO_SUPPORTED)
  set_property(TARGET MarsZ5 PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

enable_testing()
add_executable(MarsChemistryTest tests/MarsChemistryTest.cpp src/MarsChemistry.cpp)
target_include_directories(MarsChemistryTest PRIVATE include)
add_test(NAME MarsChemistryTest COMMAND MarsChemistryTest)
add_executable(MarsAtmosphereTest tests/MarsAtmosphereTest.cpp
  src/MarsAtmosphereModel.cpp src/MarsChemistry.cpp)
target_include_directories(MarsAtmosphereTest PRIVATE include)
add_test(NAME MarsAtmosphereTest COMMAND MarsAtmosphereTest)

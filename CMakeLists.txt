cmake_minimum_required(VERSION 3.5)

project(MPCNS CXX)

# C++ 
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# directories
#==================================================================
# add_subdirectory(./Z0_NULL)
option(MPCNS_BUILD_MERCURY "Configure the MercuryZ4 solver (requires PETSc)" OFF)
option(MPCNS_BUILD_MARS "Configure the MarsZ5 solver" ON)

if(MPCNS_BUILD_MERCURY)
  add_subdirectory(./Z4_Mercury)
endif()
if(MPCNS_BUILD_MARS)
  add_subdirectory(./Z5_Mars)
endif()
#==================================================================

cmake_minimum_required(VERSION 3.5)

project(MPCNS CXX)

# C++ 
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# directories
#==================================================================
add_subdirectory(./Z0_NULL)
add_subdirectory(./Z0_CoreDebug)
# add_subdirectory(./Z1_Euler)
# add_subdirectory(./Z2_MHD)
add_subdirectory(./Z4_Mercury)
# add_subdirectory(./Z5_Mercury)
#==================================================================

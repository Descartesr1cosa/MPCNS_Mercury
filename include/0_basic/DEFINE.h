#pragma once
// Control by #define

//=================================================================================================
// 1_MPCNS_Grid.cpp
// #define test_ghost_mesh

// 2_MPCNS_Output.h
#define if_Debug_output_ngg 0

// 3_field/Field_Array
#define if_Debug_Field_Array 0
//=================================================================================================

//=================================================================================================
// 1_Euler_inviscid_Roe.cpp
#define Reconstruction_Scheme 1

// 2_Euler_inviscid_Flux.cpp
#define Flux_WENO_Scheme 11
// 51 LF (52 LF-Characteristic *, only 51 & 71 are supported)
// 71 LF

// 1_NS_Viscous.h
#define Vis_Dimension 3
// 2- 2Dimension; 3- 3Dimension
//=================================================================================================

//=================================================================================================
#define Output_if_nonDimension 1
// 1-nonDimension   0-dimension

//=================================================================================================

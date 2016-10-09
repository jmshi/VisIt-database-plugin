#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file defs.hpp.in
//  \brief Template file for defs.hpp.  When the configure.py script is run, a new
//  defs.hpp file will be created (overwriting the last) from this template.  This new
//  file contains Athena++ specific cpp macros and definitions set by configure.
//======================================================================================

//--------------------------------------------------------------------------------------
// macros which define physics and algorithms

// problem generator
#define PROBLEM_GENERATOR "hgb"

// coordinate system
#define COORDINATE_SYSTEM "cartesian"

// enable shearing box? default=0 (false)
#define SHEARING_BOX 1

// non-barotropic equation of state (i.e. P not simply a func of rho)? default=1 (true)
#define NON_BAROTROPIC_EOS 0

// include magnetic fields? default=0 (false)
#define MAGNETIC_FIELDS_ENABLED 1

// include radiative transfer? default=0 (false)
#define RADIATION_ENABLED 0

// enable special or general relativity? default=0 (false)
#define RELATIVISTIC_DYNAMICS 0

// enable general relativity? default=0 (false)
#define GENERAL_RELATIVITY 0

// Riemann solver
#define RIEMANN_SOLVER "hlld"

// spatial reconstruction algorithm
#define RECONSTRUCTION_METHOD "plm"

// hydro time-integration algorithm
#define HYDRO_TIME_INTEGRATOR "vl2"

// hydro time-integration algorithm
#define COMPILED_WITH "mpicxx"

// compiler options
#define COMPILED_WITH_OPTIONS " -O3 -xhost -ipo -inline-forceinline -diag-disable 3180   -lhdf5"

// MPI parallelization (MPI_PARALLEL or NOT_MPI_PARALLEL)
//#define MPI_PARALLEL

// openMP parallelization (OPENMP_PARALLEL or NOT_OPENMP_PARALLEL)
#define NOT_OPENMP_PARALLEL

// HDF5 output (HDF5OUTPUT or NO_HDF5OUTPUT)
#define HDF5OUTPUT

//--------------------------------------------------------------------------------------
// macros associated with numerical algorithm (rarely modified)

#define NHYDRO 5
#define NFIELD 3
#define NWAVE 6
#define NGHOST 2
#define MAX_NSTEP 4

//--------------------------------------------------------------------------------------
// general purpose macros (never modified)

#define PI 3.1415926535897932
#define SQRT2 1.4142135623730951
#define ONE_OVER_SQRT2 0.7071067811865475
#define ONE_3RD 0.3333333333333333
#define TWO_3RD 0.6666666666666667
#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+36
#define SQR(x) ( (x)*(x) )
#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 )

#endif
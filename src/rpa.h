/*****************************************************************************
*
* RPA calculation for the fermionic Hubbard model
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __RPA__
#define __RPA__

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <armadillo>

typedef std::complex<double> cx_double;
typedef arma::vec vec;
typedef arma::cx_vec cx_vec;
typedef arma::mat mat;
typedef arma::cx_mat cx_mat;

// /* Using boost::filesystem */
// #include <boost/filesystem.hpp>
// typedef boost::filesystem::path path;
// typedef boost::filesystem::ofstream ofstream;
typedef std::filesystem::path path;
typedef std::ofstream ofstream;

#include "parameters.h"
#include "hoppings.h"

constexpr int NSUBL = 2;  // Number of sublattices

#endif // __RPA__

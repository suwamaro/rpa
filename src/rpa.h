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
#include <boost/filesystem.hpp>
//#include <boost/lexical_cast.hpp>
#include <armadillo>

typedef std::complex<double> cx_double;
typedef arma::vec vec;
typedef arma::cx_vec cx_vec;
typedef arma::mat mat;
typedef arma::cx_mat cx_mat;

#include "parameters.h"
#include "hoppings.h"

constexpr int NSUBL = 2;  // Number of sublattices

#endif // __RPA__

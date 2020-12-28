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

/* Uncomment if you use Boost Filesystem. */
// #define _USE_BOOST_FILESYSTEM_

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

#ifdef _USE_BOOST_FILESYSTEM_
/* Using boost::filesystem */
#include <boost/filesystem.hpp>
typedef boost::filesystem::path path;
typedef boost::filesystem::ofstream ofstream;
using boost::filesystem::exists;
#else
/* Using std::filesystem */
typedef std::filesystem::path path;
typedef std::ofstream ofstream;
using std::filesystem::exists;
#endif

using namespace std::complex_literals;

constexpr int up_spin = 1;
constexpr int down_spin = -1;

constexpr int NSUBL = 2;  // Number of sublattices

constexpr double kB = 8.617333262145e-5; // = 1380649. / 16021766340;

#endif // __RPA__

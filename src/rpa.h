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
typedef arma::sp_mat sp_mat;

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

typedef std::ostream ostream;

using namespace std::complex_literals;

constexpr int up_spin = 1;
constexpr int down_spin = -1;

constexpr int NSUBL = 2;  // Number of sublattices

constexpr double kB = 8.617333262145e-5; // = 1380649. / 16021766340;
constexpr double J_over_eV = 6.24150907446076 * 1e+18;
constexpr double e_over_hbar = 1.519267447 * 1e+15 / J_over_eV;  // (A/eV)
constexpr double planck_h = 4.135667696 * 1e-15;  // (eV/Hz)
constexpr double c_light = 2.99792458 * 1e+8;  // (m/s)
constexpr double eV_to_inv_nm = 1e-9 / (planck_h * c_light);

#endif // __RPA__

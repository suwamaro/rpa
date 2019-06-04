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
#include <boost/lexical_cast.hpp>

typedef std::complex<double> cx_double;

class hoppings {
public:
  hoppings(double theta, double phi, double t3);
  double ek1(double kx, double ky, double kz) const;
  double ek2(double kx, double ky, double kz) const;
  double ek3(double kx, double ky, double kz) const;
  double t_max() const;
  /* Parameters of the t2g three-band model */
  /* Reference: J.-M. Carter, et al., PRB 87 014433 (2013) 
                S. Mohapatra, et al., PRB 95, 094435 (2017). */
  double t1 = 0.2;
  double t2 = -0.032;
  double t3 = -0.03;     // Modified
  // double t3 = -0.02  // in Reference
  double t4 = 0.188;
  double t5 = -0.054;
  double t6 = 0;
  double tm = -0.03;
  double tmp = 0.022;
  double t1z = 0.03;
  double t4z = -0.16;
  double tmz = 0.072;
  double t6z = -0.04;
  double t2z = 0;

  /* Parameters of an effective one-band model */
  double t;
  double t_bar;
  double tp;
  double tpp;
  double tz;
  double tz_bar;
  double tzp;
};

#endif // __RPA__

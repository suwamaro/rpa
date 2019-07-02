/*****************************************************************************
*
* Class of hoppings
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"

hoppings::hoppings(double theta, double phi, double t_third){
  double ct = cos(theta);
  double ct2 = ct * ct;
  double st = sin(theta);
  double st2 = st * st;
  double cp = cos(phi);
  double sp = sin(phi);
  
  t = (t1 * st2 + (t4+t5)*ct2/2) * cp + tm * ct2 * sp;
  t_bar = (t1 * st2 - (t4+t5)*ct2/2) * sp + tm * ct2 * cp;
  tp = t2 * st2 + t6 * ct2;
  t3 = t_third;
  tpp = t3 * st2;
  tz = (t1z * st2 + t4z * ct2) * cp + tmz * ct2 * sp;
  tz_bar = (t1z * st2 - t4z * ct2) * sp + tmz * ct2 * cp;
  tzp = t2z * st2 + t6z * ct2 / 2;
  
  // // for check
  // t = 1.0;
  // t_bar = 0;
  // tp = 0;
  // tpp = 0;
  // tz = 0;
  // tz_bar = 0;
  // tzp = 0;
  
  // for check
  std::cerr << t << "  " << t_bar << "  " << tp << "  " << tpp << "  " << tz << "  " << tz_bar << "  " << tzp << std::endl;
}

double hoppings::ek1(double kx, double ky, double kz) const {
  return - 2. * ( t * cos(kx) + t * cos(ky) + 0.5 * tz * cos(kz) );
}

double hoppings::ek2(double kx, double ky, double kz) const {
  return - 2. * ( t_bar * cos(kx) + t_bar * cos(ky) + 0.5 * tz_bar * cos(kz) );
}

double hoppings::ek3(double kx, double ky, double kz) const {
  return - 2. * tp * (cos(kx+ky) + cos(kx-ky)) - 2. * tpp * (cos(2*kx) + cos(2*ky)) - 2. * tzp * cos(kz) * (cos(kx) + cos(ky));
}

double hoppings::t_max() const {
  double tmax = 0;
  tmax = std::max( tmax, sqrt(t*t + t_bar*t_bar) );
  tmax = std::max( tmax, std::abs(tp) );
  tmax = std::max( tmax, std::abs(tpp) );
  tmax = std::max( tmax, sqrt(tz*tz + tz_bar*tz_bar) );
  tmax = std::max( tmax, std::abs(tzp) );
  return tmax;
}

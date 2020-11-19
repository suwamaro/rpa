/*****************************************************************************
*
* Calculation of the single particle energy for the fermionic Hubbard model.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_single_particle_energy.h"
#include "rpa_util.h"

std::tuple<double, double> calc_single_particle_energy(hoppings const& ts, double kx, double ky, double kz, double delta){
  double ek1 = ts.ek1(kx, ky, kz);
  double ek2 = ts.ek2(kx, ky, kz);
  double ek3 = ts.ek3(kx, ky, kz);
  double Em = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double Ep = eigenenergy_HF_plus(ek1, ek2, ek3, delta);
  return std::make_tuple(Em, Ep);
}

std::tuple<double, double> calc_single_particle_energy2(hoppings2 const& ts, double kx, double ky, double kz, double delta){
  cx_double ek1 = ts.ek1(kx, ky, kz);
  cx_double ek23 = ts.ek23(kx, ky, kz);
  cx_double ekz = ts.ekz(kx, ky, kz);
  cx_double tz = ts.tz;
  double Em = eigenenergy_HF(-1, ek1, ek23, ekz, tz, kz, delta);
  double Ep = eigenenergy_HF(1, ek1, ek23, ekz, tz, kz, delta);
  return std::make_tuple(Em, Ep);
}

/*****************************************************************************
*
* Calculation of the single particle energy for the fermionic Hubbard model.
*
* Copyright (C) 2019 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __CALC_SINGLE_PARTICLE_ENERGY__
#define __CALC_SINGLE_PARTICLE_ENERGY__

#include "hoppings.h"

std::tuple<double, double> calc_single_particle_energy(hoppings const& ts, double kx, double ky, double kz, double delta);
std::tuple<double, double> calc_single_particle_energy2(hoppings2 const& ts, double kx, double ky, double kz, double delta);
  
#endif // __CALC_SINGLE_PARTICLE_ENERGY__

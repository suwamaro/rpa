/*****************************************************************************
*
* Functions for calculating the Raman scattering cross section.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _CALC_RAMAN_H
#define _CALC_RAMAN_H

#include "rpa.h"

cx_double calc_coef_eff_Raman(int L, hoppings_bilayer2 const& ts, double delta, double kx, double ky, double kz, int sigma, std::vector<BondDelta> const& bonds, BondDelta mu, BondDelta nu, vec3 const& ki, vec3 const& kf);

void calc_Raman_bilayer(path& base_dir, rpa::parameters const& pr);
void calc_Raman_bilayer2(path& base_dir, rpa::parameters const& pr);

#endif // _CALC_RAMAN_H

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

cx_double velocity_U1(hoppings2 const& ts, cx_mat const& Udg, cx_mat const& U, cx_mat const& U_bar, double kx, double ky, double kz, vec3 const& photon_q, BondDelta const& e_mu, std::vector<BondDelta> const& bonds, int sign_m, int sign_n, int sigma_m, int sigma_n);
cx_double calc_coef_eff_Raman_resonant(int L, hoppings2 const& ts, double delta, double kx, double ky, double kz, int sigma, std::vector<BondDelta> const& bonds, BondDelta mu, BondDelta nu, vec3 const& ki, vec3 const& kf, std::vector<vec3> const& occ, std::vector<vec3> const& emp);
void calc_coef_eff_Raman_nonresonant(hoppings2 const& ts, double kx, double ky, double kz, std::vector<BondDelta> const& bonds, BondDelta mu, BondDelta nu, cx_double *N);

void calc_Raman_bilayer(path& base_dir, rpa::parameters const& pr);
void calc_coef_eff_Raman_real_space(path& base_dir, rpa::parameters const& pr);

#endif // _CALC_RAMAN_H

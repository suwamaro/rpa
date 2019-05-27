/*****************************************************************************
*
* RPA calculation for the fermionic Hubbard model
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"

double energy_free_electron(double t, double mu, double kx, double ky, double kz);
double energy_free_electron(double t, double mu, double kx, double ky);
// double energy_free_electron_bilayer(hoppings const& ts, double mu, double kx, double ky, double kz);
double energy_free_electron_bilayer1(hoppings const& ts, double mu, double kx, double ky, double kz);
double energy_free_electron_bilayer2(hoppings const& ts, double mu, double kx, double ky, double kz);
double eigenenergy_HF_in(double e_free, double delta);
double eigenenergy_HF_out(double e_free, double delta);
double eigenenergy_HF_in(double ek1, double ek2, double ek3, double delta);
double eigenenergy_HF_out(double ek1, double ek2, double ek3, double delta);
double eigenenergy_HF(double e_free, double delta, double kx, double ky);
double fermi_density(double x, double kT, double mu);

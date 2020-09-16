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
double eigenenergy_HF_minus(double e_free, double delta);
double eigenenergy_HF_plus(double e_free, double delta);
double eigenenergy_HF_minus(double ek1, double ek2, double ek3, double delta);
double eigenenergy_HF_plus(double ek1, double ek2, double ek3, double delta);
cx_double bk(int spin, cx_double ek1, cx_double tz, double kz);
double zk(int spin, cx_double ek1, cx_double tz, double kz, double delta);
double zk(cx_double ek1, cx_double tz, double kz, double delta);
cx_double xk(int spin, cx_double ek1, cx_double tz, double kz, double delta);
double eigenenergy_HF(double sign, cx_double ek1, cx_double ek23, cx_double ekz, cx_double tz, double kz, double delta);
double fermi_density(double x, double kT, double mu);

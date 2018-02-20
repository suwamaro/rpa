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
double eigenenergy_HF_in(double e_free, double delta);
double eigenenergy_HF_out(double e_free, double delta);
double eigenenergy_HF(double e_free, double delta, double kx, double ky);
double fermi_density(double x, double kT, double mu);

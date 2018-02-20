/*****************************************************************************
*
* Functions for the RPA calculation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"

double energy_free_electron(double t, double mu, double kx, double ky, double kz){
  return - 2. * t * ( cos( kx ) + cos( ky ) + cos( kz ) ) - mu;
}

double energy_free_electron(double t, double mu, double kx, double ky){
  return - 2. * t * ( cos( kx ) + cos( ky ) ) - mu;
}

double eigenenergy_HF_in(double e_free, double delta){
  return - sqrt( delta * delta + e_free * e_free );
}

double eigenenergy_HF_out(double e_free, double delta){
  return sqrt( delta * delta + e_free * e_free );
}

double fermi_density(double x, double kT, double mu) {
  double alpha = (x-mu) / std::abs(kT);
  if (kT < 1e-15 || std::abs(alpha) > 20) {
    return (x < mu) ? 1.0 : 0.0;
  } else {
    return 1.0/(exp(alpha)+1.0);
  }
}

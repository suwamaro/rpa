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

double calc_tk1(hoppings const& ts, double kx, double ky, double kz){
  return ts.t * cos(kx) + ts.t * cos(ky) + 0.5 * ts.tz * cos(kz);
}

double calc_tk2(hoppings const& ts, double kx, double ky, double kz){
  return ts.t_bar * cos(kx) + ts.t_bar * cos(ky) + 0.5 * ts.tz_bar * cos(kz);  
}

double calc_tk3(hoppings const& ts, double kx, double ky, double kz){
  return ts.tp * cos(kx+ky) + ts.tp * cos(kx-ky)
    + ts.tpp * cos(2*kx) + ts.tpp * cos(2*ky)
    + ts.tzp * cos(kz) * (cos(kx) + cos(ky));
}

double energy_free_electron_bilayer1(hoppings const& ts, double mu, double kx, double ky, double kz){
  double tk1 = calc_tk1( ts, kx, ky, kz );
  double tk2 = calc_tk2( ts, kx, ky, kz );
  double tk3 = calc_tk3( ts, kx, ky, kz );
  double tk = tk3 + sqrt(tk1*tk1 + tk2*tk2);
  return - 2. * tk - mu;
}

double energy_free_electron_bilayer2(hoppings const& ts, double mu, double kx, double ky, double kz){
  double tk1 = calc_tk1( ts, kx, ky, kz );
  double tk2 = calc_tk2( ts, kx, ky, kz );
  double tk3 = calc_tk3( ts, kx, ky, kz );
  double tk = tk3 - sqrt(tk1*tk1 + tk2*tk2);
  return - 2. * tk - mu;
}

double eigenenergy_HF_minus(double e_free, double delta){
  return - sqrt( delta * delta + e_free * e_free );
}

double eigenenergy_HF_plus(double e_free, double delta){
  return sqrt( delta * delta + e_free * e_free );
}

double eigenenergy_HF_minus(double ek1, double ek2, double ek3, double delta){
  return ek3 - sqrt(delta*delta + ek1*ek1 + ek2*ek2);
}

double eigenenergy_HF_plus(double ek1, double ek2, double ek3, double delta){
  return ek3 + sqrt(delta*delta + ek1*ek1 + ek2*ek2);
}

cx_double bk(int spin, cx_double ek1, cx_double tz, double kz){
  cx_double bk0 = ek1 - tz * cos(kz);
  if (spin == up_spin) {
    return bk0;
  } else {
    return std::conj(bk0);
  }
}

double bk(cx_double ek1, cx_double tz, double kz){
  return std::abs(bk(up_spin, ek1, tz, kz));  // Spin does not matter.
}

double zk(int spin, cx_double ek1, cx_double tz, double kz, double delta){
  return delta / sqrt( delta*delta + std::norm(bk(spin,ek1,tz,kz)));
}

double zk(cx_double ek1, cx_double tz, double kz, double delta){
  return zk(up_spin, ek1, tz, kz, delta);  // Spin does not matter.
}

double zk_over_delta(int spin, cx_double ek1, cx_double tz, double kz, double delta){
  return 1. / sqrt( delta*delta + std::norm(bk(spin,ek1,tz,kz)));
}

double zk_over_delta(cx_double ek1, cx_double tz, double kz, double delta){
  return zk_over_delta(up_spin, ek1, tz, kz, delta);
}

cx_double xk(int spin, cx_double ek1, cx_double tz, double kz, double delta){
  cx_double bki = bk(spin, ek1, tz, kz);
  double bk_abs = std::abs(bki);
  if ( bk_abs < 1e-12 ) {
    return 1.0;
  } else {
    return bki / std::abs(bki);
  }
}

double eigenenergy_HF(double sign, cx_double ek1, cx_double ek23, cx_double ekz, cx_double tz, double kz, double delta){
  return std::real(ek23) + std::real(ekz) + sign * sqrt(delta*delta + std::norm(bk(up_spin,ek1,tz,kz)));
}

double fermi_energy(double x, double kT, double mu) {
  double alpha = (x-mu)/std::abs(kT);
  if (kT < 1e-15 || std::abs(alpha) > 40) {
    return (x < mu) ? (x-mu) : 0.0;
  }
  else {
    if ( alpha >= 0 ) {
      return - kT * log(1. + exp(-alpha));
    } else {
      return x - mu - kT * log(1. + exp(alpha));
    }
  }
}

double fermi_density(double x, double kT, double mu) {    
  double alpha = (x-mu) / std::abs(kT);
  if (kT < 1e-15 || std::abs(alpha) > 40.) {
    if (std::abs(x-mu) < 1e-12) {
      return 0.5;
    } else if (x < mu) {
      return 1.0;
    } else {
      return 0.0;
    }    
  } else {
    if ( alpha > 0 ) {
      return exp(-alpha) / (exp(-alpha) + 1.);
    } else {
      return 1. / (exp(alpha) + 1.);
    }    
  }
}

double compressibility(double x, double kT, double mu) {    
  double alpha = (x-mu) / std::abs(kT);
  if (kT < 1e-15 || std::abs(alpha) > 40.) {
    return 0; /* Not considering the delta peak */
  } else {
    if ( alpha > 0 ) {
      return exp(-alpha) / std::pow(exp(-alpha) + 1., 2) / kT;
    } else {
      return exp(alpha) / std::pow(exp(alpha) + 1., 2) / kT;
    }    
  }
}

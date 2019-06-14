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

// double energy_free_electron_bilayer(hoppings const& ts, double mu, double kx, double ky, double kz){
//   // 0.5 for z direction



//   double sign = - 1.;
//   if ( std::abs(kz) < 1e-10 ) {
//     // if ( std::abs(kx) + std::abs(ky) <= M_PI + 1e-10 ) {
//     if ( std::abs(kx) + std::abs(ky) < M_PI ) {
      
//       sign = 1.;
//       // sign = - 1.;
      
//     } else {
//       sign = - 1.;
//       // sign = 1.;
      
//     }
//   } else if ( std::abs(kz + M_PI) < 1e-10 ) {
//     // if ( std::abs(kx) + std::abs(ky) < M_PI ) {
//     if ( std::abs(kx) + std::abs(ky) <= M_PI + 1e-10 ) {
      
//       sign = 1.;
//       // sign = - 1.;
      
//     } else {
//       sign = - 1.;
//       // sign = 1.;
      
//     }
//   } else { /* Not being considered for the moment. */ }
//   double tk = tk3 + sign * sqrt(tk1*tk1 + tk2*tk2);  
//   // return - 2. * tk - mu;
//   return - 2. * tk1 - mu;
// }

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

double fermi_density(double x, double kT, double mu) {
  double alpha = (x-mu) / std::abs(kT);
  if (kT < 1e-15 || std::abs(alpha) > 20) {
    return (x < mu) ? 1.0 : 0.0;
  } else {
    return 1.0/(exp(alpha)+1.0);
  }
}

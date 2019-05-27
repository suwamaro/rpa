/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_gap.h"
#include "rpa_util.h"

double denominator_in(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_in(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = delta + Ek - ek3;
  double denom = sqrt( g*g + std::norm(nume) );
  return denom;
}

cx_double calc_ak_in(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_in(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = delta + Ek - ek3;
  double denom = sqrt( g*g + std::norm(nume) );

  // // for check
  // std::cerr << "nume = " << nume << "  g = " << g << "  denom = " << denom << std::endl;
  
  return nume / denom;  
}

cx_double calc_ak_out(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_out(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = delta + Ek - ek3;
  double denom = sqrt( g*g + std::norm(nume) );
  return nume / denom;  
}

double calc_bk_in(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_in(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = delta + Ek - ek3;
  double denom = sqrt( g*g + std::norm(nume) );
  return g / denom;  
}

double calc_bk_out(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_out(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = delta + Ek - ek3;
  double denom = sqrt( g*g + std::norm(nume) );
  return g / denom;  
}


double calc_bk_up_in(double e_free, double delta){
  double Ek = eigenenergy_HF_in( e_free, delta );
  return  e_free / sqrt( ( Ek - delta ) * ( Ek - delta ) + e_free * e_free );
}

double calc_bk_up_out(double e_free, double delta){
  double Ek = eigenenergy_HF_out( e_free, delta );
  
    /* Special treatment on the Fermi surface */
  double e_eps = 1e-12;
  if ( std::abs( e_free ) <= e_eps ) { return 1.; }
  else {
    return  e_free / sqrt( ( Ek - delta ) * ( Ek - delta ) + e_free * e_free );
  }
}

double calc_ak_up_in(double e_free, double delta){
  double bk = calc_bk_up_in( e_free, delta );
  return sqrt( 1. - bk * bk );
}

double calc_ak_up_out(double e_free, double delta){
  double bk = calc_bk_up_out( e_free, delta );
  return sqrt( 1. - bk * bk );
}

double calc_bk_down_in(double e_free, double delta){
  return calc_ak_up_in( e_free, delta );
}

double calc_bk_down_out(double e_free, double delta){
  return calc_ak_up_out( e_free, delta );
}

double calc_ak_down_in(double e_free, double delta){
  return calc_bk_up_in( e_free, delta );
}

double calc_ak_down_out(double e_free, double delta){
  return calc_bk_up_out( e_free, delta );
}

cx_double larger_eigenvalue(cx_double A, cx_double B, cx_double D){
  return 0.5 * ( A + D + sqrt( std::conj( A - D ) * ( A - D ) + 4. * std::conj(B) * B ) );
}

double wave_vector_in_BZ(double k){
  /* Returning - M_PI <= k < M_PI */
  double k_eps = 1e-12;
  while( k > M_PI - k_eps ) k -= 2. * M_PI;
  while( k < - M_PI - k_eps ) k += 2. * M_PI;
  if ( std::abs( k ) < k_eps ) k = 0;
  if ( std::abs( k + M_PI ) < k_eps ) k = - M_PI;
  return k;
}

void add_to_sus_mat(cx_double& A, cx_double& B, cx_double& D, double e_free, double e_free2, double delta, cx_double omega){
  double e_eps = 1e-12;
  double E1 = eigenenergy_HF_out( e_free2, delta );
  double E2 = eigenenergy_HF_in( e_free, delta );
  cx_double diff_E1 = E1 - E2 + omega;
  // cx_double diff_E2 = E1 - E2 - std::conj(omega);
  cx_double diff_E2 = E1 - E2 - omega;
  
  double ak_up_in = calc_ak_up_in( e_free, delta );
  double ak_q_up_out = calc_ak_up_out( e_free2, delta );
  double ak_down_in = calc_ak_down_in( e_free, delta );
  double ak_q_down_out = calc_ak_down_out( e_free2, delta );
  double bk_up_in = calc_bk_up_in( e_free, delta );
  double bk_q_up_out = calc_bk_up_out( e_free2, delta );
  double bk_down_in = calc_bk_down_in( e_free, delta );
  double bk_q_down_out = calc_bk_down_out( e_free2, delta );
  
  double factor = 1.;
  
  /* Taking into account a half of the contribution from the Fermi surface */
  if ( std::abs( e_free ) <= e_eps ) {
    factor = 0.5;
  }
  
  A += factor * ( ak_up_in * ak_up_in * ak_q_down_out * ak_q_down_out / diff_E1 + ak_down_in * ak_down_in * ak_q_up_out * ak_q_up_out / diff_E2 );
  B += factor * ( ak_up_in * bk_up_in * ak_q_down_out * bk_q_down_out / diff_E1 + ak_down_in * bk_down_in * ak_q_up_out * bk_q_up_out / diff_E2 );
  D += factor * ( bk_up_in * bk_up_in * bk_q_down_out * bk_q_down_out / diff_E1 + bk_down_in * bk_down_in * bk_q_up_out * bk_q_up_out / diff_E2 );
}

// void add_to_sus_mat(cx_double& A, cx_double& B, cx_double& D, double e_free, double e_free2, double delta, cx_double omega){
//   double e_eps = 1e-12;
//   double E1 = eigenenergy_HF_out( e_free2, delta );
//   double E2 = eigenenergy_HF_in( e_free, delta );
//   cx_double diff_E1 = E1 - E2 + omega;
//   // cx_double diff_E2 = E1 - E2 - std::conj(omega);
//   cx_double diff_E2 = E1 - E2 - omega;
  
//   double ak_up_in = calc_ak_up_in( e_free, delta );
//   double ak_q_up_out = calc_ak_up_out( e_free2, delta );
//   double ak_down_in = calc_ak_down_in( e_free, delta );
//   double ak_q_down_out = calc_ak_down_out( e_free2, delta );
//   double bk_up_in = calc_bk_up_in( e_free, delta );
//   double bk_q_up_out = calc_bk_up_out( e_free2, delta );
//   double bk_down_in = calc_bk_down_in( e_free, delta );
//   double bk_q_down_out = calc_bk_down_out( e_free2, delta );
  
//   double factor = 1.;
  
//   /* Taking into account a half of the contribution from the Fermi surface */
//   if ( std::abs( e_free ) <= e_eps ) {
//     factor = 0.5;
//   }
  
//   A += factor * ( ak_up_in * ak_up_in * ak_q_down_out * ak_q_down_out / diff_E1 + ak_down_in * ak_down_in * ak_q_up_out * ak_q_up_out / diff_E2 );
//   B += factor * ( ak_up_in * bk_up_in * ak_q_down_out * bk_q_down_out / diff_E1 + ak_down_in * bk_down_in * ak_q_up_out * bk_q_up_out / diff_E2 );
//   D += factor * ( bk_up_in * bk_up_in * bk_q_down_out * bk_q_down_out / diff_E1 + bk_down_in * bk_down_in * bk_q_up_out * bk_q_up_out / diff_E2 );
// }

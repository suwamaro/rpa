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

/* Reference: A. Singh and Z. Tesanovic, PRB 41, 11457 (1990) */

bool diagonal_H(double ek1, double ek2) {  
  return ek1 * ek1 + ek2 * ek2 < 1e-12;
}
double diagonal(double ek3, double delta, double E){
  return E - ek3 - delta;
}
cx_double off_diagonal(double ek1, double ek2){
  cx_double comp(ek1, - ek2);
  return comp;  
}
cx_double calc_bk_up_in_minus(double ek1, double ek2, double ek3, double delta){
  double E = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double diag = diagonal(ek3, delta, E);
  cx_double off_diag = off_diagonal(ek1, ek2);
  double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
  return off_diag / denom;
}
cx_double calc_ak_up_in_minus(double ek1, double ek2, double ek3, double delta){
  double E = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double diag = diagonal(ek3, delta, E);
  cx_double off_diag = off_diagonal(ek1, ek2);
  double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
  return diag / denom;  
}
cx_double calc_ak_down_in_minus(double ek1, double ek2, double ek3, double delta){
  return calc_bk_up_in_minus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_down_in_minus(double ek1, double ek2, double ek3, double delta){
  return calc_ak_up_in_minus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_up_in_plus(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_bk_down_in_minus( - ek1, - ek2, ek3, delta ));  /* k + pi */
}
cx_double calc_ak_up_in_plus(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_ak_down_in_minus( - ek1, - ek2, ek3, delta ));  /* k + pi */  
}
cx_double calc_ak_down_in_plus(double ek1, double ek2, double ek3, double delta){
  return calc_bk_up_in_plus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_down_in_plus(double ek1, double ek2, double ek3, double delta){
  return calc_ak_up_in_plus(ek1, ek2, ek3, delta);
}

double calc_bk_up_in(double e_free, double delta){
  double Ek = eigenenergy_HF_minus( e_free, delta );
  return  e_free / sqrt( ( Ek - delta ) * ( Ek - delta ) + e_free * e_free );
}

double calc_bk_up_out(double e_free, double delta){
  double Ek = eigenenergy_HF_plus( e_free, delta );
  
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
  double E1 = eigenenergy_HF_plus( e_free2, delta );
  double E2 = eigenenergy_HF_minus( e_free, delta );
  cx_double diff_E1 = E1 - E2 + omega;
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

void add_to_sus_mat2(hoppings const& ts, double mu, cx_double& A, cx_double& B, cx_double& C, cx_double& D, double qx, double qy, double qz, double kx, double ky, double kz, double delta, cx_double omega, bool zz){
  double eps = 1e-12;
  
  double ek1 = ts.ek1(kx, ky, kz);
  double ek2 = ts.ek2(kx, ky, kz);
  double ek3 = ts.ek3(kx, ky, kz);
  double Ek = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  
  /* Checking if the eigenenergy is below the chemical potential. */
  double mu_free = 0;  /* Assume at half filling */
  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
  if ( e_free > mu_free + eps ) return;
    
  /* Prefactor */
  double factor = 1.;
  if ( zz ) { factor = 0.25; }
  
  /* On the zone boundary */
  if ( std::abs(e_free - mu_free) < eps ) {
    factor *= 0.5;
  }
  
  double diff_x = wave_vector_in_BZ( kx - qx );
  double diff_y = wave_vector_in_BZ( ky - qy );
  double diff_z = wave_vector_in_BZ( kz - qz );
  
  double ek_q1 = ts.ek1(diff_x, diff_y, diff_z);
  double ek_q2 = ts.ek2(diff_x, diff_y, diff_z);
  double ek_q3 = ts.ek3(diff_x, diff_y, diff_z);
  double Ek_q = eigenenergy_HF_plus(ek_q1, ek_q2, ek_q3, delta);
  
  cx_double diff_E1 = Ek_q - Ek + omega;
  cx_double diff_E2 = Ek_q - Ek - omega;
  
  cx_double ak_up = calc_ak_up_in_minus(ek1, ek2, ek3, delta);
  cx_double ak_down = calc_ak_down_in_minus(ek1, ek2, ek3, delta);
  cx_double bk_up = calc_bk_up_in_minus(ek1, ek2, ek3, delta);
  cx_double bk_down = calc_bk_down_in_minus(ek1, ek2, ek3, delta);

  /* For k-q inside the magnetic BZ, the negative sign for b comes from k -> k + pi. */
  /* For k-q outside the magnetic BZ, the negative sign for b comes from the Fourier transform. */
  double sign_A = 1.;
  double sign_B = - 1.;
  
  cx_double ak_q_up = sign_A * calc_ak_up_in_plus( ek_q1, ek_q2, ek_q3, delta);
  cx_double ak_q_down = sign_A * calc_ak_down_in_plus( ek_q1, ek_q2, ek_q3, delta);  
  cx_double bk_q_up = sign_B * calc_bk_up_in_plus( ek_q1, ek_q2, ek_q3, delta);
  cx_double bk_q_down = sign_B * calc_bk_down_in_plus( ek_q1, ek_q2, ek_q3, delta);
  
  if ( zz ) {
    A += factor * ( std::norm(ak_up) * std::norm(ak_q_up) / diff_E1 + std::norm(ak_up) * std::norm(ak_q_up) / diff_E2 );
    A += factor * ( std::norm(ak_down) * std::norm(ak_q_down) / diff_E1 + std::norm(ak_down) * std::norm(ak_q_down) / diff_E2 );    
    
    B += factor * ( std::conj(ak_up) * bk_up * ak_q_up * std::conj(bk_q_up) / diff_E1 + std::conj(ak_up) * bk_up * ak_q_up * std::conj(bk_q_up) / diff_E2 );
    B += factor * ( std::conj(ak_down) * bk_down * ak_q_down * std::conj(bk_q_down) / diff_E1 + std::conj(ak_down) * bk_down * ak_q_down * std::conj(bk_q_down) / diff_E2 );    
    
    C += factor * ( ak_up * std::conj(bk_up) * std::conj(ak_q_up) * bk_q_up / diff_E1 + ak_up * std::conj(bk_up) * std::conj(ak_q_up) * bk_q_up / diff_E2 );
    C += factor * ( ak_down * std::conj(bk_down) * std::conj(ak_q_down) * bk_q_down / diff_E1 + ak_down * std::conj(bk_down) * std::conj(ak_q_down) * bk_q_down / diff_E2 );    

    D += factor * ( std::norm(bk_up) * std::norm(bk_q_up) / diff_E1 + std::norm(bk_up) * std::norm(bk_q_up) / diff_E2 );
    D += factor * ( std::norm(bk_down) * std::norm(bk_q_down) / diff_E1 + std::norm(bk_down) * std::norm(bk_q_down) / diff_E2 );    
  } else {
    A += factor * ( std::norm(ak_up) * std::norm(ak_q_down) / diff_E1 + std::norm(ak_down) * std::norm(ak_q_up) / diff_E2 );
    
    B += factor * ( std::conj(ak_up) * bk_up * ak_q_down * std::conj(bk_q_down) / diff_E1 + std::conj(ak_down) * bk_down * ak_q_up * std::conj(bk_q_up) / diff_E2 );
    
    C += factor * ( ak_up * std::conj(bk_up) * std::conj(ak_q_down) * bk_q_down / diff_E1 + ak_down * std::conj(bk_down) * std::conj(ak_q_up) * bk_q_up / diff_E2 );

    D += factor * ( std::norm(bk_up) * std::norm(bk_q_down) / diff_E1 + std::norm(bk_down) * std::norm(bk_q_up) / diff_E2 );
  }
}

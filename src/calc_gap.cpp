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

double denominator_in(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_in(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = Ek - delta - ek3;
  double denom = sqrt( g*g + std::norm(nume) );
  return denom;
}

double denominator_out(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_out(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = Ek - delta - ek3;
  double denom = sqrt( g*g + std::norm(nume) );
  return denom;
}

cx_double calc_bk_up_in(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_in(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = Ek - delta - ek3;
  double denom = sqrt( g*g + std::norm(nume) );

  // // for check
  // std::cerr << Ek << " " << ek1 << " " << g << std::endl;
  
  return nume / denom;
}

cx_double calc_bk_up_out(double ek1, double ek2, double ek3, double delta){
  double Ek = eigenenergy_HF_out(ek1, ek2, ek3, delta);
  cx_double nume(ek1, ek2);
  double g = Ek - delta - ek3;
  double denom = sqrt( g*g + std::norm(nume) );
  /* Special case */
  if ( std::abs(denom) < 1e-12 ) {
    return 1.;
  } else {  
    return nume / denom;
  }
}

cx_double calc_ak_up_in(double ek1, double ek2, double ek3, double delta){
  return sqrt( 1 - std::norm( calc_bk_up_in(ek1, ek2, ek3, delta) ) );
  
}

cx_double calc_ak_up_out(double ek1, double ek2, double ek3, double delta){
  return sqrt( 1 - std::norm( calc_bk_up_out(ek1, ek2, ek3, delta) ) );
}

cx_double calc_bk_down_in(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_ak_up_in(ek1, ek2, ek3, delta));
}

cx_double calc_bk_down_out(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_ak_up_out(ek1, ek2, ek3, delta));
}

cx_double calc_ak_down_in(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_bk_up_in(ek1, ek2, ek3, delta));
}

cx_double calc_ak_down_out(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_bk_up_out(ek1, ek2, ek3, delta));
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

  //   // for check
  // std::cerr << "e_free = " << e_free << "  e_free2 = " << e_free2 << std::endl;
  // std::cerr << std::real(bk_up_in) << " " << std::real(bk_q_up_out) << std::endl;
  // std::cerr << std::real(ak_up_in) << " " << std::real(ak_q_down_out) << " " << std::real(ak_down_in) << " " << std::real(ak_q_up_out) << " " << std::real(bk_up_in) << " " << std::real(bk_q_down_out) << " " << std::real(bk_down_in) << " " << std::real(bk_q_up_out) << std::endl;
  A += factor * ( ak_up_in * ak_up_in * ak_q_down_out * ak_q_down_out / diff_E1 + ak_down_in * ak_down_in * ak_q_up_out * ak_q_up_out / diff_E2 );
  B += factor * ( ak_up_in * bk_up_in * ak_q_down_out * bk_q_down_out / diff_E1 + ak_down_in * bk_down_in * ak_q_up_out * bk_q_up_out / diff_E2 );
  D += factor * ( bk_up_in * bk_up_in * bk_q_down_out * bk_q_down_out / diff_E1 + bk_down_in * bk_down_in * bk_q_up_out * bk_q_up_out / diff_E2 );
}

void add_to_sus_mat2(hoppings const& ts, double mu, cx_double& A, cx_double& B, cx_double& D, double qx, double qy, double qz, double kx, double ky, double kz, double delta, cx_double omega){
  double e_eps = 1e-12;
  
  /* Checking if k is inside/outside the BZ. */
  // double k_len = std::abs(kx) + std::abs(ky) + std::abs(kz);  
  double k_len = std::abs(kx) + std::abs(ky);

  /* Outside the BZ */
  if ( k_len - M_PI >=  e_eps ) return;

  double diff_x = wave_vector_in_BZ( kx - qx );
  double diff_y = wave_vector_in_BZ( ky - qy );
  double diff_z = wave_vector_in_BZ( kz - qz );
  
  double ek1 = ts.ek1(kx, ky, kz);
  double ek2 = ts.ek2(kx, ky, kz);
  double ek3 = ts.ek3(kx, ky, kz);
  
  double ek_q1 = ts.ek1(diff_x, diff_y, diff_z);
  double ek_q2 = ts.ek2(diff_x, diff_y, diff_z);
  double ek_q3 = ts.ek3(diff_x, diff_y, diff_z);
  
  double Ek = eigenenergy_HF_in(ek1, ek2, ek3, delta);
  double Ek_q = eigenenergy_HF_out(ek_q1, ek_q2, ek_q3, delta);

  // // for check
  // std::cerr << qx << " " << qy << " " << qz << " " << kx << " " << ky << " " << kz << "  " << Ek << "  " << Ek_q;
  // if ( std::abs(Ek-mu) < e_eps ) std::cerr << "  0\n";
  // else if ( Ek - mu > e_eps ) std::cerr << "  positive\n";
  // else std::cerr << "  negative\n";
  
  /* Checking if the eigenenergy is below the chemical potential. */
  if ( Ek > mu + e_eps ) return;

  // // for check
  // std::cerr << "Considering..."  << std::endl;
  
  /* Prefactor */
  double factor = 1.;
  
  /* On the zone boundary */
  if ( std::abs(k_len - M_PI) < e_eps ) {
    factor = 0.5;
  }
		
  // /* Checking if the denominators are zero. */
  // double denom_in = denominator_in(ek1, ek2, ek3, delta);
  // double denom_out = denominator_out(ek_q1, ek_q2, ek_q3, delta);  

  cx_double ak_up_in = calc_ak_up_in(ek1, ek2, ek3, delta);
  cx_double ak_q_up_out = calc_ak_up_out(ek_q1, ek_q2, ek_q3, delta);
  cx_double ak_down_in = calc_ak_down_in(ek1, ek2, ek3, delta);
  cx_double ak_q_down_out = calc_ak_down_out(ek_q1, ek_q2, ek_q3, delta);
  cx_double bk_up_in = calc_bk_up_in(ek1, ek2, ek3, delta);
  cx_double bk_q_up_out = calc_bk_up_out(ek_q1, ek_q2, ek_q3, delta);
  cx_double bk_down_in = calc_bk_down_in(ek1, ek2, ek3, delta);
  cx_double bk_q_down_out = calc_bk_down_out(ek_q1, ek_q2, ek_q3, delta);    
  
  // /* Special case */
  // if ( std::abs(denom_in) < e_eps ) {
  //   ak = 1.;
  //   bk = 0.;
  // } else {
    // ak = calc_ak_in(ek1, ek2, ek3, delta);
    // bk = calc_bk_in(ek1, ek2, ek3, delta);    
  // }
  
  // if ( std::abs(denom_out) < e_eps ) {
  //   ak_q = 1.;
  //   bk_q = 0.;
  // } else {
    // ak_q = calc_ak_out(ek_q1, ek_q2, ek_q3, delta);    
    // bk_q = calc_bk_out(ek_q1, ek_q2, ek_q3, delta);
  // }
    
  cx_double diff_E1 = Ek_q - Ek + omega;
  // cx_double diff_E2 = Ek_q - Ek - std::conj(omega);
  cx_double diff_E2 = Ek_q - Ek - omega;

  // // for check
  // std::cerr << qx << " " << qy << " " << kx << " " << ky << " ";
  // std::cerr << std::real(bk_up_in) << " " << std::real(bk_q_up_out) << std::endl;  
  // std::cerr << std::real(ak_up_in) << " " << std::real(ak_q_down_out) << " " << std::real(ak_down_in) << " " << std::real(ak_q_up_out) << " " << std::real(bk_up_in) << " " << std::real(bk_q_down_out) << " " << std::real(bk_down_in) << " " << std::real(bk_q_up_out) << std::endl;
    // std::cerr << qx << " " << qy << " " << kx << " " << ky << " " << diff_x << " " << diff_y << " " << std::real(ak_up_in) << " " << std::real(ak_q_down_out) << " " << std::real(ak_down_in) << " " << std::real(ak_q_up_out) << " " << std::real(bk_up_in) << " " << std::real(bk_q_down_out) << " " << std::real(bk_down_in) << " " << std::real(bk_q_up_out) << std::endl;  
  // std::cerr << qx << " " << qy << " " << qz << " " << kx << " " << ky << " " << kz << " " << diff_x << " " << diff_y << " " << diff_z << " " << std::real(ak_up_in) << " " << std::real(ak_q_down_out) << " " << std::real(ak_down_in) << " " << std::real(ak_q_up_out) << " " << std::real(bk_up_in) << " " << std::real(bk_q_down_out) << " " << std::real(bk_down_in) << " " << std::real(bk_q_up_out) << std::endl;

  // std::cerr << "delta = " << delta << std::endl;
  // std::cerr << "eks" << std::endl;
  // std::cerr << ek1 << " " << ek2 << " " << ek3 << " " << ek_q1 << " " << ek_q2 << " " << ek_q3 << std::endl;
  // // std::cerr << "denom_in = " << denom_in << std::endl;
  // // std::cerr << "denom_out = " << denom_out << std::endl;
  // std::cerr << Ek << " " << Ek_q << " " << ak << " " << ak_q << " " << bk << " " << bk_q << std::endl;
  
  // A += factor * ( std::norm(ak) * bk_q * bk_q / diff_E1 + bk * bk * std::norm(ak_q) / diff_E2 );
  // // B += factor * ak * bk * std::conj(ak_q) * bk_q * ( 1. / diff_E1 + 1. / diff_E2 );
  // B += factor * ak * bk * ak_q * bk_q * ( 1. / diff_E1 + 1. / diff_E2 );
  
  // D += factor * ( bk * bk * std::norm(ak_q) / diff_E1 + std::norm(ak) * bk_q * bk_q / diff_E2 );

  A += factor * ( std::norm(ak_up_in) * std::norm(ak_q_down_out) / diff_E1 + std::norm(ak_down_in) * std::norm(ak_q_up_out) / diff_E2 );
  B += factor * ( ak_up_in * bk_up_in * ak_q_down_out * bk_q_down_out / diff_E1 + ak_down_in * bk_down_in * ak_q_up_out * bk_q_up_out / diff_E2 );
  D += factor * ( std::norm(bk_up_in) * std::norm(bk_q_down_out) / diff_E1 + std::norm(bk_down_in) * std::norm(bk_q_up_out) / diff_E2 );  
}

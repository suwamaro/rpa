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
#include "BinarySearch.h"

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

double larger_eigenvalue(double A, double B, double D){
  return 0.5 * ( A + D + sqrt( ( A - D ) * ( A - D ) + 4. * B * B ) );
}

double wave_vector_in_BZ(double k){
  double k_eps = 1e-12;
  while( k > M_PI - k_eps ) k -= 2. * M_PI;
  while( k < - M_PI - k_eps ) k += 2. * M_PI;
  if ( std::abs( k ) < k_eps ) k = 0;
  if ( std::abs( k + M_PI ) < k_eps ) k = - M_PI;
  return k;
}

double calc_eigval(int L, double t, double mu, double U, double delta, double qx, double qy, double omega){
  double k1 = 2. * M_PI / (double)L;
  
  double A = 0, B = 0, D = 0;
  for(int x=-L/2; x < L/2; x++){
    double kx = k1 * x;
    for(int y=-L/2; y < L/2; y++){
      double ky = k1 * y;

      double e_free = energy_free_electron( t, mu, kx, ky );
      double e_eps = 1e-12;
      if ( e_free < e_eps ) {
	/* Summing up over all k inside the Brillouin zone. */
	double diff_x = wave_vector_in_BZ( kx - qx );
	double diff_y = wave_vector_in_BZ( ky - qy );	
	double e_free2 = energy_free_electron( t, mu, diff_x, diff_y );
	double E1 = eigenenergy_HF_out( e_free2, delta );
	double E2 = eigenenergy_HF_in( e_free, delta );
	double diff_E1 = E1 - E2 + omega;
	double diff_E2 = E1 - E2 - omega;

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
    }
  }

  int n_sites = L * L;
  A *= 2. / (double)n_sites;
  B *= 2. / (double)n_sites;
  D *= 2. / (double)n_sites;
  return larger_eigenvalue( A, B, D );
}

double calc_gap(int L, double t, double mu, double U, double delta, double qx, double qy){
  double target = 1. / U;
  double omega = 0.01;
  
  using std::placeholders::_1;
  auto calc_ev = std::bind( calc_eigval, L, t, mu, U, delta, qx, qy, _1 );

  BinarySearch bs;
  double omega_delta = 0.001;
  // bs.find_solution( omega, target, calc_ev );
  // bs.find_solution( omega, target, calc_ev, true, omega_delta, 0, 1., true );
  bs.find_solution( omega, target, calc_ev, true, omega_delta );
  
  return omega;
}

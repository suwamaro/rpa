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

double calc_eigval_cubic(int L, double t, double mu, double U, double delta, double qx, double qy, double qz, double omega){
  double k1 = 2. * M_PI / (double)L;
  
  cx_double A = 0, B = 0, D = 0;
  for(int x=-L/2; x < L/2; x++){
    double kx = k1 * x;
    for(int y=-L/2; y < L/2; y++){
      double ky = k1 * y;

      for(int z=-L/2; z < L/2; z++){
	double kz = k1 * z;
	
	double e_free = energy_free_electron( t, mu, kx, ky, kz );
	double e_eps = 1e-12;
	if ( e_free < e_eps ) {
	  /* Summing up over all k inside the Brillouin zone. */
	  double diff_x = wave_vector_in_BZ( kx - qx );
	  double diff_y = wave_vector_in_BZ( ky - qy );
	  double diff_z = wave_vector_in_BZ( kz - qz );	
	  double e_free2 = energy_free_electron( t, mu, diff_x, diff_y, diff_z );
	  add_to_sus_mat( A, B, D, e_free, e_free2, delta, omega );
	}
      }
    }
  }

  int n_sites = L * L * L;
  A *= 2. / (double)n_sites;
  B *= 2. / (double)n_sites;
  D *= 2. / (double)n_sites;
  return std::real( larger_eigenvalue( A, B, D ) );
}

double calc_gap_cubic(int L, double t, double mu, double U, double delta, double qx, double qy, double qz){
  double target = 1. / U;
  double omega = 0.01;
  
  using std::placeholders::_1;
  auto calc_ev = std::bind( calc_eigval_cubic, L, t, mu, U, delta, qx, qy, qz, _1 );

  BinarySearch bs;
  double omega_delta = 0.001;
  // bs.find_solution( omega, target, calc_ev );
  // bs.find_solution( omega, target, calc_ev, true, omega_delta, 0, 1., true );
  bs.find_solution( omega, target, calc_ev, true, omega_delta );
  
  return omega;
}

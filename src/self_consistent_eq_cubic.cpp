/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "self_consistent_eq.h"
#include "rpa_util.h"
#include "BinarySearch.h"

double self_consistent_eq_cubic(int L, double t, double mu, double delta){
  /* Monotonically decreasing as a function of delta */
  double k1 = 2. * M_PI / (double)L;
  double sum = 0;
  for(int x=-L/2; x < L/2; x++){
    double kx = k1 * x;
    for(int y=-L/2; y < L/2; y++){
      double ky = k1 * y;

      // for check
      for(int z=-L/2; z < L/2; z+=L/2){	      
      // for(int z=-L/2; z < L/2; z++){
	double kz = k1 * z;
      
	double e_free = energy_free_electron( t, mu, kx, ky, kz );
	double e_eps = 1e-12;

	/* Summing up over all k inside the Brillouin zone. */
	if ( e_free < e_eps ) {
	  double factor = 1.;
	  
	  /* Taking into account a half of the contribution from the Fermi surface */
	  if ( std::abs( e_free ) <= e_eps ) {
	    factor = 0.5;
	  }
	  
	  sum += factor / eigenenergy_HF_out( e_free, delta );
	}
      }
    }
  }

  // for check
  int n_sites = L * L * 2;
  // int n_sites = L * L * L;  
  return sum / (double)( n_sites );
}

double solve_self_consistent_eq_cubic(int L, double t, double mu, double U){
  std::cout << "Finding a self-consistent solution for U=" << U << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_cubic, L, t, mu, _1 );

  BinarySearch bs;
  bs.find_solution( delta, target, scc );
  return delta;
}


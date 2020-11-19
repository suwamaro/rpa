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
#include "calc_gap.h"
#include "BinarySearch.h"

double self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double delta){
  /* Monotonically decreasing as a function of delta */
  double k1 = 2. * M_PI / (double)L;
  double sum = 0;
  for(int x=-L/2; x < L/2; x++){
    double kx = k1 * x;
    for(int y=-L/2; y < L/2; y++){
      double ky = k1 * y;
      for(int z=-L/2; z < L/2; z++){	
	double kz = k1 * z;
      
	double e_free = energy_free_electron( ts.t, mu, kx, ky, kz );
	double e_eps = 1e-12;

	/* Summing up over all k inside the Brillouin zone. */
	if ( e_free < e_eps ) {
	  double factor = 1.;

	  /* Taking into account a half of the contribution from the Fermi surface */
	  if ( std::abs( e_free ) <= e_eps ) {
	    factor = 0.5;
	  }
	  
	  /* Sum of the Fourier transformed hoppings */
	  double ek1 = ts.ek1(kx, ky, kz);
	  double ek2 = ts.ek2(kx, ky, kz);
	  double ek3 = ts.ek3(kx, ky, kz);
	  
	  cx_double ak_up = calc_ak_up_in_minus(ek1, ek2, ek3, delta);
	  cx_double ak_down = calc_ak_down_in_minus(ek1, ek2, ek3, delta);
	  
	  sum += factor * (std::norm(ak_up) - std::norm(ak_down));	  
	  // sum += factor / eigenenergy_HF_plus( e_free, delta );
	}
      }
    }
  }

  int n_sites = L * L * L;
  sum /= (double)( n_sites * delta );

  /* The self-consistent condition: sum = 1/U */
  return sum;
}

double solve_self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double U){
  std::cout << "Finding a self-consistent solution for U=" << U << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_cubic, L, ts, mu, _1 );

  BinarySearch bs;
  bool sol_found = bs.find_solution( delta, target, scc );
  if ( sol_found ) {
    return delta;
  } else {
    return 0;
  }    
}


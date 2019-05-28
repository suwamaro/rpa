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

double self_consistent_eq_bilayer(int L, hoppings const& ts, double mu, double delta){
  /* Monotonically decreasing as a function of delta */
  double k1 = 2. * M_PI / (double)L;
  double sum = 0;

  // for check
  for(int z=0; z < 1; z++){
  // for(int z=-1; z < 1; z++){
    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;      
	double e_eps = 1e-12;
	
	// // for check
	// std::cerr << "kx = " << kx << "  ky = " << ky << "  kz = " << kz << std::endl;

	/* Checking if k is inside/outside the BZ. */
	double k_len = std::abs(kx) + std::abs(ky);

	/* Outside the BZ */
	if ( k_len - M_PI >=  e_eps ) continue;

	/* Prefactor */
	double factor = 1.;
	
	/* On the zone boundary */
	if ( std::abs(k_len - M_PI) < e_eps ) {
	  factor = 0.5;
	}

	/* Sum of the Fourier transformed hoppings */
	double ek1 = ts.ek1(kx, ky, kz);
	double ek2 = ts.ek2(kx, ky, kz);
	double ek3 = ts.ek3(kx, ky, kz);

	// // for check
	// std::cerr << "ek1 = " << ek1 << "  ek2 = " << ek2 << "  ek3 = " << ek3 << std::endl;
	
	// /* Checking if the denominator is zero. */
	// double denom_in = denominator_in(ek1, ek2, ek3, delta);

	// std::cerr << "denom_in = " << denom_in << std::endl;

	// /* Special case */
	// if ( std::abs(denom_in) < e_eps ) {
	//   sum += factor;
	//   continue;
	// }

	/* Eigenenergy inside the BZ */
	double Ek_in = eigenenergy_HF_in(ek1, ek2, ek3, delta) - mu;

	/* Summing up for Ek < EF. */
	if ( Ek_in < e_eps ) {
	  cx_double ak_up = calc_ak_up_in(ek1, ek2, ek3, delta);
	  cx_double ak_down = calc_ak_down_in(ek1, ek2, ek3, delta);

	  // // for check
	  // std::cerr << "ak = " << ak << "   bk = " << bk << std::endl;

	  /* |a|^2 - |b|^2 */
	  sum += factor * (std::norm(ak_up) - std::norm(ak_down));
	}
      } /* end for y */
    } /* end for x */
  } /* end for z */


  // for check
  int n_sites = L * L;
  // int n_sites = L * L * 2;
  
  sum /= (double)( n_sites * delta );

  // // for check
  // std::cerr << delta << "  " << sum << std::endl;

  /* The self-consistent condition: sum = 1/U */
  return sum;
}

double solve_self_consistent_eq_bilayer(int L, hoppings const& ts, double mu, double U){
  std::cout << "Finding a self-consistent solution for U=" << U << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_bilayer, L, ts, mu, _1 );

  BinarySearch bs;
  bs.find_solution( delta, target, scc );
  return delta;
}


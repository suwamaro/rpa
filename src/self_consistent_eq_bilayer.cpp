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

/* Creating an instance */
SelfConsistentIntegrandBilayer scib;

int sc_integrand_bilayer(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return scib.calc(ndim, xx, ncomp, ff, userdata);
}

double self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double delta){
  double eps = 1e-12;
  
  /* Monotonically decreasing as a function of delta */
  double k1 = 2. * M_PI / (double)L;
  double sum = 0;

  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;      
	double e_eps = 1e-12;
	
	/* Checking if k is inside/outside the BZ. */
	double mu_free = 0;  /* Assume at half filling */
	double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
	if ( e_free > mu_free + eps ) continue;

	/* Prefactor */
	double factor = 1.;
	
	/* On the zone boundary */
	if ( std::abs(e_free - mu_free) < eps ) {	
	  factor = 0.5;
	}

	/* Sum of the Fourier transformed hoppings */
	double ek1 = ts.ek1(kx, ky, kz);
	double ek2 = ts.ek2(kx, ky, kz);
	double ek3 = ts.ek3(kx, ky, kz);

	cx_double ak_up = calc_ak_up_in_minus(ek1, ek2, ek3, delta);
	cx_double ak_down = calc_ak_down_in_minus(ek1, ek2, ek3, delta);

	sum += factor * (std::norm(ak_up) - std::norm(ak_down));
      } /* end for y */
    } /* end for x */
  } /* end for z */

  int n_sites = L * L * 2;  
  sum /= (double)( n_sites * delta );

  /* The self-consistent condition: sum = 1/U */
  return sum;
}

double solve_self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double U){
  std::cout << "Finding a self-consistent solution for U=" << U << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_bilayer, L, ts, _1 );

  BinarySearch bs;
  bs.find_solution( delta, target, scc );
  return delta;
}

double self_consistent_eq_bilayer2(int L, hoppings_bilayer2& ts, double delta, CubaParam const& cbp, bool continuous_k){  
  /* Monotonically decreasing as a function of delta */
  double sum = 0;
  
  if ( continuous_k ) {
    /* Changing the parameters */
    scib.set_parameters(&ts, delta);

    /* For Cuba */
    double epsrel = 1e-8;
    double epsabs = 1e-10;    
    int ncomp = 1;
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, sc_integrand_bilayer, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);    
    // for(int comp = 0; comp < cbp.NCOMP; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);

    sum = integral[0] / ( 4. * delta );
  } else {
    double eps = 1e-12;
    double k1 = 2. * M_PI / (double)L;
  
    for(int z=0; z < 2; z++){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; x++){
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; y++){
	  double ky = k1 * y;      
	
	  /* Checking if k is inside/outside the BZ. */
	  double mu_free = 0;  /* Assume at half filling */
	  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
	  if ( e_free > mu_free + eps ) continue;

	  /* Prefactor */
	  double factor = 1.;
	
	  /* On the zone boundary */
	  if ( std::abs(e_free - mu_free) < eps ) {	
	    factor = 0.5;
	  }

	  /* Sum of the Fourier transformed hoppings */
	  cx_double ek1 = ts.ek1(kx, ky, kz);
	  double zki = zk(ek1, ts.tz, kz, delta);
	  
	  sum += factor * zki;
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    sum /= (double)( n_sites * delta );    
  }  

  /* The self-consistent condition: sum = 1/U */
  return sum;
}

double solve_self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double U, CubaParam const& cbp, bool continuous_k){
  
  std::cout << "Finding a self-consistent solution for U=" << U << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_bilayer2, L, ts, _1, cbp, continuous_k );

  BinarySearch bs;
  bs.find_solution( delta, target, scc );
  return delta;
}


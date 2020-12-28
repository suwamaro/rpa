/*****************************************************************************
*
* Functions for calculating the critical T.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "find_critical_T.h"
#include "BinarySearch.h"
#include "calc_chemical_potential.h"

/* Instantiation */
FindTcIntegrandBilayer ftib;

/* Integrand */
int find_critical_T_bilayer_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return ftib.calc(ndim, xx, ncomp, ff, userdata);
}

double self_consistent_eq_T_bilayer2(int L, hoppings_bilayer2 const& ts, double filling, double T, CubaParam const& cbp, bool continuous_k){
  double sum = 0;

  /* Changing the parameters */
  double mu = calc_chemical_potential_bilayer3(L, ts, filling, T, 0, cbp, continuous_k, false);  
  ftib.set_parameters(ts, T, mu);
    
  if ( continuous_k ) {  
    /* For Cuba */
    double epsrel = 1e-10;
    double epsabs = 1e-10;
    int ncomp = 1;    
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, find_critical_T_bilayer_integrand, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    // for(int comp = 0; comp < ncomp; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);      

    sum = integral[0] / 4.;
  } else { /* Integral for a finite-size system */
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

	  double qvec[3] = {kx, ky, kz};
	  sum += factor * ftib.integrand(qvec);
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    sum /= (double)(n_sites);    
  }

  return sum;
}

double find_critical_T_bilayer(rpa::parameters const& pr){
  /* Getting parameters */  
  int L = pr.L;
  double U = pr.U;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  std::cout << "Finding a self-consistent solution for U = " << U << std::endl;
  double target = 1. / U;
  double T = 1e-12;
    
  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_T_bilayer2, L, std::ref(*ts), filling, _1, std::ref(cbp), continuous_k );

  BinarySearch bs(continuous_k);
  bs.set_x_MIN(0);
  double T_delta = 100. * U / kB;
  bool verbose = false;
  bool sol_found = bs.find_solution( T, target, scc, true, T_delta, verbose );

  double Tc = 0;
  if ( sol_found ) {
    Tc = T;
  } else {
    std::cerr << "Not found." << std::endl;
    Tc = 0;
  }
  
  return Tc;
}

void find_critical_T_bilayer_output(path& base_dir, rpa::parameters const& pr){
  std::cout << "Finding the critical T..." << std::endl;
  double Tc = find_critical_T_bilayer(pr);
  
  /* Precision */
  int prec = 12;
  
  /* Output */
  std::cout << "Tc = " << Tc << std::endl;
  
  ofstream out_Tc;
  out_Tc.open( base_dir / "critical-T.text");
  out_Tc << "Tc = " << std::setprecision(prec) << std::setw( prec ) << Tc << std::endl;
  out_Tc.close();  
}

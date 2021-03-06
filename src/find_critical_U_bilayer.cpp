/*****************************************************************************
*
* Functions for calculating the critical U.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "find_critical_U.h"

/* Instantiation */
FindUcIntegrandBilayer fuib;

/* Integrand */
int find_critical_U_bilayer_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return fuib.calc(ndim, xx, ncomp, ff, userdata);
}

double find_critical_U_bilayer(rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Setting the parameters */
  fuib.set_parameters(*ts, filling);  
  
  double sum = 0;  
  if ( continuous_k ) {
    /* Parameters for Cuba */
    CubaParam cbp(pr);
  
    /* For Cuba */
    double epsrel = 1e-10;
    double epsabs = 1e-10;
    int ncomp = 1;    
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, find_critical_U_bilayer_integrand, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
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
	  sum += factor * fuib.integrand(qvec);	  
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    sum /= (double)(n_sites);    
  }

  /* Uc */  
  double Uc = 1. / sum;
  return Uc;
}

void find_critical_U_bilayer_output(path& base_dir, rpa::parameters const& pr){
  std::cout << "Finding the critical U..." << std::endl;
  double Uc = find_critical_U_bilayer(pr);
  
  /* Precision */
  int prec = 12;
  
  /* Output */
  std::cout << "Uc = " << Uc << std::endl;
  
  ofstream out_Uc;
  out_Uc.open( base_dir / "critical-U.text");
  out_Uc << "Uc = " << std::setprecision(prec) << std::setw( prec ) << Uc << std::endl;
  out_Uc.close();  
}

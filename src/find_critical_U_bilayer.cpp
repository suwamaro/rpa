/*****************************************************************************
*
* Functions for calculating the spectrum
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

/* Member functions of FindUcIntegrandBilayer */
void FindUcIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h){
  hb_ = h;
}

int FindUcIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
  /* Reset */
  ff[0] = 0;
  
  /* Wavenumbers */
  double k1 = xx[0] * 2 * M_PI;
  double k2 = xx[1] * 2 * M_PI;
  
  double kx = 0.5 * (k2 + k1);
  double ky = 0.5 * (k2 - k1);
  
  /* Sum over kz */
  for(int z=0; z < 2; z++){       
    double kz = M_PI * z;	  
    
    cx_double ek1 = ts()->ek1(kx, ky, kz);
    ff[0] += 1. / std::abs(bk(up(), ek1, ts()->tz, kz));  // Spin does not matter.
  }
  
  return 0;   
}

void find_critical_U_bilayer(path& base_dir, rpa::parameters const& pr){
  std::cout << "Finding the critical U..." << std::endl;
  
  /* Getting parameters */
  int L = pr.L;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  double t1 = pr.t1;
  double t2 = pr.t2;
  double t3 = pr.t3;
  double tz1 = pr.t4;
  double tz2 = pr.t5;
  double phase1 = pr.phase1;
  double phasez = pr.phase4;
  
  using namespace std::complex_literals;  
  // cx_double t1_cx( pr.t1, pr.t1_bar );
  // cx_double tz_cx( pr.t4, pr.t4_bar );
  cx_double t1_cx = t1*exp(1i*phase1*0.5);
  cx_double tz1_cx = tz1*exp(1i*phasez*0.5);
  
  /* Hopping class */
  std::unique_ptr<hoppings_bilayer2> ts;
  ts = hoppings_bilayer2::mk_bilayer2(t1_cx, t2, t3, tz1_cx, tz2);
    
  double sum = 0;  
  if ( continuous_k ) {
    /* Changing the parameters */
    fuib.set_parameters(*ts);

    /* Parameters for Cuba */
    CubaParam cbp(pr);
  
    /* For Cuba */
    double epsrel = 1e-8;
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

	  cx_double ek1 = ts->ek1(kx, ky, kz);
	  cx_double bki = bk(up(), ek1, ts->tz, kz);  // Spin does not matter.
	  sum += factor / std::abs(bki);
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    sum /= (double)(n_sites);    
  }

  /* Uc */  
  double Uc = 1. / sum;
  
  /* Precision */
  int prec = 12;
  
  /* Output */
  std::cout << "Uc = " << Uc << std::endl;
  
  ofstream out_Uc;
  out_Uc.open( base_dir / "critical-U.text");
  out_Uc << "Uc = " << std::setw( prec ) << Uc << std::endl;
  out_Uc.close();  
}

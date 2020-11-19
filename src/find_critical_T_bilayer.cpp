/*****************************************************************************
*
* Functions for calculating the critical T.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "calc_single_particle_energy.h"
#include "find_critical_T.h"
#include "BinarySearch.h"

/* Instantiation */
FindTcIntegrandBilayer ftib;

/* Integrand */
int find_critical_T_bilayer_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return ftib.calc(ndim, xx, ncomp, ff, userdata);
}

/* Member functions of FindTcIntegrandBilayer */
void FindTcIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double T){
  hb_ = h;
  set_T(T);
}

int FindTcIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
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

    /* Fermi density */
    double n_minus = 1.0;
    double n_plus = 0.0;
    if ( kB * T() < 1e-15 ) {
      n_minus = 1.0;
      n_plus = 0.0;
    } else {
      double Em, Ep;
      std::tie(Em, Ep) = calc_single_particle_energy2(*ts(), kx, ky, kz, 0);  // delta == 0
      double mu = 0.5 * ( Em + Ep );
      n_minus = fermi_density(Em, kB*T(), mu);
      n_plus = fermi_density(Ep, kB*T(), mu);	    
    }
	  
    ff[0] += (n_minus - n_plus) / std::abs(bk(up_spin, ek1, ts()->tz, kz));  // Spin does not matter.
  }
  
  return 0;   
}

double self_consistent_eq_T_bilayer2(int L, hoppings_bilayer2 const& ts, double T, CubaParam const& cbp, bool continuous_k){
  double sum = 0;  
  if ( continuous_k ) {
    /* Changing the parameters */
    ftib.set_parameters(ts, T);
  
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

	  cx_double ek1 = ts.ek1(kx, ky, kz);
	  cx_double bki = bk(up_spin, ek1, ts.tz, kz);  // Spin does not matter.

	  /* Fermi density */
	  double n_minus = 1.0;
	  double n_plus = 0.0;
	  if ( kB * T < 1e-15 ) {
	    n_minus = 1.0;
	    n_plus = 0.0;
	  } else {
	    double Em, Ep;
	    std::tie(Em, Ep) = calc_single_particle_energy2(ts, kx, ky, kz, 0);  // delta == 0
	    double mu = 0.5 * ( Em + Ep );
	    n_minus = fermi_density(Em, kB*T, mu);
	    n_plus = fermi_density(Ep, kB*T, mu);	    
	  }
	  
	  sum += factor * (n_minus - n_plus) / std::abs(bki);
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    sum /= (double)(n_sites);    
  }

  return sum;
}

void find_critical_T_bilayer(path& base_dir, rpa::parameters const& pr){
  std::cout << "Finding the critical T..." << std::endl;
  
  /* Getting parameters */  
  int L = pr.L;
  double U = pr.U;
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

  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  std::cout << "Finding a self-consistent solution for U = " << U << std::endl;
  double target = 1. / U;
  double T = 0.1 * U;
    
  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_T_bilayer2, L, std::ref(*ts), _1, std::ref(cbp), continuous_k );

  BinarySearch bs;
  bs.set_x_MIN(0);
  bool sol_found = bs.find_solution( T, target, scc );

  double Tc = 0;
  if ( sol_found ) {
    Tc = T;
  }
  
  /* Precision */
  int prec = 12;
  
  /* Output */
  std::cout << "Tc = " << Tc << std::endl;
  
  ofstream out_Tc;
  out_Tc.open( base_dir / "critical-T.text");
  out_Tc << "Tc = " << std::setprecision(prec) << std::setw( prec ) << Tc << std::endl;
  out_Tc.close();  
}

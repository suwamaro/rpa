/*****************************************************************************
*
* Functions for calculating electronic current.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <complex>
#ifdef WITH_OpenMP
#include <omp.h>
#endif
#include "calc_intensity.h"
#include "calc_spectrum.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"
#include "calc_chemical_potential.h"
#include "calc_current.h"
#include "BinarySearch.h"

/* Creating an instance */
CurrentIntegrandBilayer cuib;

/* Member functions of CurrentIntegrandBilayer */
cx_double CurrentIntegrandBilayer::integrand(double kx, double ky, double kz){    
  double kvec[3] = {kx, ky, kz};
  set_current_matrix(kvec, dir());

  cx_double current = 0;
  for(std::size_t n=0; n < 4; n++){
    double E = mfb_.E(n);
    double f = fermi_density(E, kT(), mu());
    if ( f > 1e-14 ) {
      current += arma::cdot(mfb_.U()->col(n), J_ * mfb_.U()->col(n)) * f;
    }
  }
  
  return current;
}

void CurrentIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double kT, double mu, double delta, int dir, double *Qvec){
  hb_ = h;
  kT_ = kT;
  mu_ = mu;
  delta_ = delta;
  dir_ = dir;
  for(int d=0; d<3; d++){ Qvec_[d] = Qvec[d]; }
}


int CurrentIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  cx_double sum = 0;
  
  /* Wavenumbers */
  double k1 = xx[0] * 2 * M_PI;
  double k2 = xx[1] * 2 * M_PI;
  
  double kx = 0.5 * (k2 + k1);
  double ky = 0.5 * (k2 - k1);
  
  /* Sum over kz */
  for(int z=0; z < 2; z++){       
    double kz = M_PI * z;	  
    sum += integrand(kx, ky, kz);
  }

  ff[0] = std::real(sum);
  ff[1] = std::imag(sum);
  
  return 0;   
}

void CurrentIntegrandBilayer::set_current_matrix(double* kvec, int dir){
  mfb_.set_parameters(hb_, delta_, kvec);
  mfb_.set_eigen_energies();
  mfb_.set_eigen_vectors();

  cx_double t = 0;
  double r_diff[3] = {0, 0, 0};
  if ( dir == x_dir() ) {
    r_diff[0] = 1.;    
    t = - hb_.t;
  } else if ( dir == y_dir() ) {
    r_diff[1] = 1.;
    t = - hb_.t;
  } else if ( dir == z_dir() ) {
    r_diff[2] = 1.;
    t = - hb_.tz;
  } else {
    std::cerr << "\"dir\" has to be x, y, or z." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  double inner_prod = 0;
  for(int d=0; d < 3; d++){ inner_prod += kvec[d] * r_diff[d]; }
  cx_double pure_i(0,1);
  cx_double phase = std::exp( pure_i * inner_prod );
  cx_double elem = pure_i * (t * phase - t * std::conj(phase));
  
  double Qinner_prod = 0;
  for(int d=0; d < 3; d++){ Qinner_prod += Qvec_[d] * r_diff[d]; }
  cx_double Qphase = std::exp( pure_i * Qinner_prod );

  /* Assume the U(1) symmetry */  
  J_.zeros();
  J_(2,0) = elem;
  J_(0,2) = std::conj(elem) * Qphase;
  J_(3,1) = std::conj(elem);
  J_(1,3) = elem * Qphase;
}

int current_integrand_bilayer(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return cuib.calc(ndim, xx, ncomp, ff, userdata);
}

void calc_current_bilayer(path& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  int Lk = pr.Lk;
  double U = pr.U;
  double filling = pr.filling;
  double T = pr.T;
  bool continuous_k = pr.continuous_k;

  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);
  
  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  /* Calculate the chemical potential and the charge gap. */
  double delta = solve_self_consistent_eq_bilayer2( L, *ts, U, filling, T, cbp, continuous_k );  
  std::cout << "delta = " << delta << std::endl;  
  /* Assume that mu does not depend on L for integral over continuous k. */
  double ch_gap, mu;
  std::tie(ch_gap, mu) = calc_charge_gap_bilayer( L, *ts, delta );  /* Finite size */  
  // double mu = calc_chemical_potential_bilayer2( base_dir, L, *ts, delta );  /* Finite size */

  /* Precision */
  int prec = 15;
  
  /* Output */
  ofstream out_gap;
  out_gap.open(base_dir / "gap.text");
  out_gap << "Charge gap = " << ch_gap << std::endl;
  out_gap << "mu = " << mu << std::endl;
  out_gap.close();
    
  /* Calculating electronic current */
  cx_double current[3] = {0, 0, 0};

  /* Output */
  ofstream out_current;
  out_current.open(base_dir / "electronic_current.out");
  
  double Qvec[3];
  for(int Qidx=0; Qidx < 2; Qidx++){
    if ( Qidx == 0 ) {
      for(int d=0; d < 3; d++){ Qvec[d] = 0; }
    } else {
      for(int d=0; d < 3; d++){ Qvec[d] = M_PI; }      
    }
    
    for(int dir=0; dir < 3; dir++){    
      /* Setting the parameters */  
      cuib.set_parameters(*ts, kB*T, mu, delta, dir, Qvec);
    
      /* Calculating the weight distribution */
      if ( continuous_k ) {
	/* For Cuba */
	double epsrel = 1e-10;
	double epsabs = 1e-10;    
	int ncomp = 2;
	int nregions, neval, fail;
	cubareal integral[ncomp], error[ncomp], prob[ncomp];

	/* Cuhre */
	Cuhre(cbp.NDIM, ncomp, current_integrand_bilayer, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
	// for check
	printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);    
	// for(int comp = 0; comp < ncomp; comp++ )
	//   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
	// 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);

	current[dir] = { integral[0] / 2., integral[1] / 2. };
      } else {
	double eps = 1e-12;
	double k1 = 2. * M_PI / (double)L;

	cx_double current_k = 0;
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
	  
	      current_k += factor * cuib.integrand(kx, ky, kz);
	    } /* end for y */
	  } /* end for x */
	} /* end for z */

	int n_sites = L * L * 2;
	int n_unit_cell = n_sites / 2;
	current_k /=  (double)n_unit_cell;
	current[dir] = current_k;
      }
    }

    /* Output */
    out_current << std::setprecision(8);    
    out_current << "# q = ( " << Qvec[0] << " " << Qvec[1] << " " << Qvec[2] << " )" << std::endl;
    out_current << std::setprecision(prec);
    for(int dir=0; dir < 3; dir++){
      out_current << dir << std::setw(prec+10) << current[dir] * e_over_hbar << std::endl;
    }
    out_current << std::endl;
  }  /* End for Qidx */
  
  out_current.close();
}

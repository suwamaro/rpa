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
#include "rpa_util.h"
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
cx_double CurrentIntegrandBilayer::integrand(double *ks){
  if ( no_contribution() ) { return 0.0; }
  
  set_variables(ks);
    
  cx_double current1 = 0, current2 = 0;
  for(int sigma1: {up_spin, down_spin}){    
    for(int sigma2: {up_spin, down_spin}){
      int g1s1 = sublattice_spin_index(gamma1(), sigma1);
      int g1s2 = sublattice_spin_index(gamma1(), sigma2);
      int g2s1 = sublattice_spin_index(gamma2(), sigma1);
      int g2s2 = sublattice_spin_index(gamma2(), sigma2);

      cx_double t1 = hopping_U1(gamma2(), gamma1(), sigma1, sigma2);
      double ks1[3] = {ks[0] + Qvec_[0]/2, ks[1] + Qvec_[1]/2, ks[2] + Qvec_[2]/2};  
      double inner_prod1 = ks1[0] * bond_.x + ks1[1] * bond_.y + ks1[2] * bond_.z;
      cx_double phase1 = std::exp(1i * inner_prod1);
      current1 += 1i * t1 * UfUd_(g1s2, g2s1) * epsilon() * phase1;
      
      cx_double t2 = hopping_U1(gamma2(), gamma1(), sigma2, sigma1);
      double ks2[3] = {- ks[0] + Qvec_[0]/2, - ks[1] + Qvec_[1]/2, - ks[2] + Qvec_[2]/2};  
      double inner_prod2 = ks2[0] * bond_.x + ks2[1] * bond_.y + ks2[2] * bond_.z;
      cx_double phase2 = std::exp(1i * inner_prod2);      
      current2 += 1i * std::conj(t2) * UfUd_(g2s2, g1s1) * epsilon() * phase2;
    }
  }
  
  return current1 - current2;  
}

void CurrentIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double kT, double mu, double delta, BondDelta bd, int g1, int g2, double *Qvec){
  hb_ = h;
  kT_ = kT;
  mu_ = mu;
  delta_ = delta;
  bond_ = bd;
  gamma1_ = g1;
  gamma2_ = g2;
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
    double ks[3] = {kx, ky, kz};
    sum += integrand(ks);
  }

  ff[0] = std::real(sum);
  ff[1] = std::imag(sum);
  
  return 0;
}

bool CurrentIntegrandBilayer::no_contribution() const {
  if ( bond_ == BondDelta(1, 0, 0) || bond_ == BondDelta(0, 1, 0) ) {
    if ( gamma1() == gamma2() ) {
      return true;
    }
  } else if ( bond_ == BondDelta(1, 1, 0) || bond_ == BondDelta(1, -1, 0) ) {
    if ( gamma1() != gamma2() ) {
      return true;
    }
  } else if ( bond_ == BondDelta(0, 0, 1) ) {
    if ( gamma1() == gamma2() ) {
      return true;
    }    
  } else {
    std::cerr << "This type of \"BondDelta\" is not supported in " << __func__ << "." << std::endl;
    std::exit(EXIT_FAILURE);
  }  
  return false;
}

void CurrentIntegrandBilayer::set_variables(double* ks){
  /* Setting the unitary matrix for k. */
  mfb_.set_parameters(hb_, delta(), ks);
  mfb_.set_eigen_energies();
  mfb_.set_eigen_vectors();

  /* Setting UfUd matrix. */
  UfUd_.zeros(4, 4);
  for(int i=0; i < 4; ++i){
    for(int j=0; j < 4; ++j){
      for(int n=0; n < 4; ++n){
	double E = mfb_.E(n);
	double f = fermi_density(E, kT(), mu());
	UfUd_(i,j) += (*mfb_.U())(i,n) * f * std::conj((*mfb_.U())(j,n));
      }
    }
  }

  /* Setting the epsilon. */
  if ( Qvec_[0] == 0 && Qvec_[1] == 0 && Qvec_[2] == 0 ) {
    epsilon_ = 1.0;
  } else {
    /* Assume Qvec == (pi, pi, pi) */
    epsilon_ = gamma2() == 0 ? - 1.0 : 1.0;
  }
}

cx_double CurrentIntegrandBilayer::hopping_U1(int g1, int g2, int sigma1, int sigma2){
  /* U(1) symmetric case */
  if ( sigma1 != sigma2 ) { return 0; }
  
  if ( bond_ == BondDelta(1, 0, 0) || bond_ == BondDelta(0, 1, 0) ) {
    if ( g1 == g2 ) {
      return 0;
    } else {
      return - hb_.t;
    }
  } else if ( bond_ == BondDelta(1, 1, 0) || bond_ == BondDelta(1, -1, 0) ) {
    if ( g1 == g2 ) {
      return - hb_.tp;
    } else {
      return 0;
    }
  } else if ( bond_ == BondDelta(0, 0, 1) ) {
    if ( g1 == g2 ) {
      return 0;
    } else {
      if ( g2 == 0 ) {
	if ( sigma2 == up_spin ) {
	  return - hb_.tz;
	} else {
	  return - std::conj(hb_.tz);
	}
      } else {
	if ( sigma2 == up_spin ) {
	  return - std::conj(hb_.tz);
	} else {
	  return - hb_.tz;
	}
      }
    }    
  } else {
    std::cerr << "This type of \"BondDelta\" is not supported in " << __func__ << "." << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
    
  /* Output */
  ofstream out_current;
  out_current.open(base_dir / "electronic_current.out");

  std::vector<BondDelta> bonds {
    {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, -1, 0},  {0, 0, 1}
  };
  
  double Qvec[3];
  for(int Qidx=0; Qidx < 2; Qidx++){
    if ( Qidx == 0 ) {
      for(int d=0; d < 3; d++){ Qvec[d] = 0; }  // Q=(0,0,0)      
    } else {
      for(int d=0; d < 3; d++){ Qvec[d] = M_PI; }  // Q=(pi,pi,pi)
    }

    out_current << std::setprecision(8);
    out_current << std::endl << "# q = ( " << Qvec[0] << " " << Qvec[1] << " " << Qvec[2] << " )" << std::endl;
    out_current << "# delta_x  delta_y  delta_z  gamma1  gamma2   Re[J]   Im[J]" << std::endl;    
    out_current << std::setprecision(prec);
  
    for(BondDelta bond: bonds){
      std::vector<std::pair<int,int>> sublattice_pairs{{0,0}, {0,1}, {1,0}, {1,1}};
      for(std::pair<int,int> gamma_pair12: sublattice_pairs){		      
	int g1 = gamma_pair12.first;
	int g2 = gamma_pair12.second;

	/* Calculating electronic current */
	cx_double current = 0;
  
	/* Setting the parameters */  
	cuib.set_parameters(*ts, kB*T, mu, delta, bond, g1, g2, Qvec);
    
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

	  current = { integral[0] / 2., integral[1] / 2. };
	} else {
	  double eps = 1e-12;
	  double k1 = 2. * M_PI / (double)L;

	  for(int z=0; z < 2; z++){    
	    double kz = M_PI * z;
	    for(int x=-L/2; x < L/2; x++){
	      double kx = k1 * x;
	      for(int y=-L/2; y < L/2; y++){
		double ky = k1 * y;

		/* Checking if the wavevector is inside the BZ. */
		double factor = BZ_factor_square_half_filling(kx, ky);
		if ( std::abs(factor) < 1e-12 ) { continue; }

		double ks[3] = {kx, ky, kz};
		current += factor * cuib.integrand(ks);		
	      } /* end for y */
	    } /* end for x */
	  } /* end for z */
	}

	/* Output */
	out_current << bond.x << "  " << bond.y << "  " << bond.z << "  " << g1 << "  " << g2 << std::setw(prec+10) << std::real(current) << std::setw(prec+10) << std::imag(current) << std::endl;
      }  /* end for gamma_pair12 */
    }  /* end for bond */
  }  /* End for Qidx */
  
  out_current.close();
}

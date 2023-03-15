/*****************************************************************************
*
* Functions for calculating the critical U.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "calc_chemical_potential.h"
#include "self_consistent_eq.h"
#include "find_critical_U.h"

extern SelfConsistentIntegrand2Bilayer sci2b;

/* Instantiation */
FindUcIntegrandBilayer fuib;

/* Integrand */
int find_critical_U_bilayer_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return fuib.calc(ndim, xx, ncomp, ff, userdata);
}

double find_critical_U_bilayer(rpa::parameters const& pr){
  std::cout << "Finding the critical U..." << std::endl;
  
  /* Getting parameters */
  int L = pr.L;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Parameters for Cuba */
  CubaParam cbp(pr);
    
  /* Calculating the chemical potential */
  double delta = 0.0;   // Zero order parameter
  double T = 0.0;   // Zero temperature
  double mu = calc_chemical_potential_bilayer3(L, *ts, filling, T, delta, cbp, continuous_k, false);
  // double mu = calc_chemical_potential_bilayer3(L, *ts, filling, T, delta, cbp, continuous_k, true);

  /* Setting the parameters */
  fuib.set_parameters(*ts, delta, mu);  
  
  double sum = 0;  
  if ( continuous_k ) {  
    /* For Cuba */
    int ncomp = 1;    
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, find_critical_U_bilayer_integrand, cbp.userdata, cbp.nvec, cbp.epsrel, cbp.epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
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

void check_mean_field_eq_bilayer(path const& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Parameters for Cuba */
  CubaParam cbp(pr);

  /* Precision */
  int prec = 12;
  int sw = prec + 10;
    
  /* Output file */
  double shift = std::pow(10, 6);  
  double t4r = std::round(pr.t4*shift)/shift;
  std::string ofn = "mean_field_eq-t4_"+std::to_string(t4r)+".text";
  ofstream out_mf(base_dir/ofn);
  out_mf << "#  delta   mu   func   energy" << std::endl;
  out_mf << std::setprecision(prec);
    
  double T = 0.0;   // Zero temperature
  
  /* Calculating the chemical potential */
  double delta_max = 1.0;
  double delta_delta = 0.01;
  for(double delta=0.; delta <= delta_max; delta += delta_delta){
    double mu = calc_chemical_potential_bilayer3(L, *ts, filling, T, delta, cbp, continuous_k, false);

    /* Set parameters */
    double U = 0.;   // Unused
    sci2b.set_parameters(pr, L, *ts, U, filling, T, delta, mu, continuous_k, true);

    double f = sci2b.calc_mean_field_function();
    double E = sci2b.calc_energy();
    out_mf << delta << std::setw(sw) << mu << std::setw(sw) << f << std::setw(sw) << E << std::endl;
  }

  /* Closing the output file. */
  out_mf.close();    
}

std::tuple<double, double, double, double, double, double> calc_total_energies(rpa::parameters const& pr, double U, double delta_i, double mu_i, bool set_initial_values){
  /* Getting parameters */
  int L = pr.L;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Parameters for Cuba */
  CubaParam cbp(pr);
    
  /* For the disordered state */
  double delta = 0;
  double T = 0.0;
  double mu = calc_chemical_potential_bilayer3(L, *ts, filling, T, delta, cbp, continuous_k, false);
  sci2b.set_parameters(pr, L, *ts, U, filling, T, delta, mu, continuous_k, true);
  double F1 = sci2b.calc_energy();

  /* For the ordered state */
  if (set_initial_values){
    delta = delta_i;
    mu = mu_i;
  } else {
    delta = 0.1 * U;
  }
  bool non_zero_delta = true;
  std::tie(delta, mu) = solve_self_consistent_eqs_bilayer2(pr, L, *ts, U, filling, T, continuous_k, non_zero_delta, delta, mu);
  // std::tie(delta, mu) = solve_self_consistent_eqs_bilayer(pr, L, *ts, U, filling, T, continuous_k, delta, mu);  
  sci2b.set_input(delta, mu);
  double F2 = sci2b.calc_energy();
  double diff = sci2b.calc_diff(delta, mu);

  /* Calculating the charge gap. */
  double ch_gap, mu0;
  std::tie(ch_gap, mu0) = calc_charge_gap_bilayer(L, *ts, delta);
  
  return std::make_tuple(F1, F2, delta, mu, ch_gap, diff);
}

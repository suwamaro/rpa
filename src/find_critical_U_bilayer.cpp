/*****************************************************************************
*
* Functions for calculating the critical U.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "BinarySearch.h"
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
  double delta_delta = 0.001;
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

std::tuple<double, double, double, double, double, double> calc_total_energies_bilayer(rpa::parameters const& pr, double U, double delta_i, double mu_i, bool set_initial_values){
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

double calc_grand_potential_bilayer(SelfConsistentIntegrand2Bilayer& sc, double delta){
  double mu = sc.calc_chemical_potential(delta);
  sc.set_input(delta, mu);
  return sc.calc_energy();
}

bool find_first_order_transition_point_bilayer(rpa::parameters const& pr, double& U, double delta_i){
  std::cout << "Finding a finite delta to produce the same grand potential with the disordered state..." << std::endl;

  /* Getting parameters */
  int L = pr.L;
  double T = pr.T;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Parameters for Cuba */
  CubaParam cbp(pr);

  /* Setting the parameters */
  double delta0 = 0.;
  double mu0 = 0.;
  bool non_zero_delta = true;  
  sci2b.set_parameters(pr, L, *ts, U, filling, T, delta0, mu0, continuous_k, non_zero_delta);
  
  /* Calculating the grand potential of the disordered state. */
  double F1 = calc_grand_potential_bilayer(sci2b, delta0);
  std::cout << "F1 = " << F1 << std::endl;
  
  /* Finding the order parameter (delta) to match the disordered grand potential. */
  double target = F1;
  using std::placeholders::_1;
  auto eq = std::bind(calc_grand_potential_bilayer, std::ref(sci2b), _1);

  BinarySearch bs(continuous_k);
  bs.set_x_MIN(pr.find_U1st_delta_min);
  bs.set_x_MAX(pr.find_U1st_delta_max);
  double x_delta = pr.find_U1st_delta_delta;
  
  std::cout << "Binary search to find a solution." << std::endl;  
  double delta = delta_i;
  bool additive = true;
  bool sol_found = bs.find_solution(delta, target, eq, additive, x_delta, pr.debug_mode);
  if (!sol_found) {
    return false;
  }

  /* Calculating the corresponding value of U. */
  double mu = sci2b.calc_chemical_potential(delta);
  sci2b.set_input(delta, mu);
  double U_inv = sci2b.calc_mean_field_function();
  sci2b.set_U(1./U_inv);
  U = 1./U_inv;
  
  /* Output */
  double F2 = sci2b.calc_energy();
  double Fdiff = F2 - F1;
  double scdiff = sci2b.calc_diff();
  std::cout << "delta = " << delta << "   mu = " << mu << "   U_inv = " << U_inv << "   U1st = " << 1. / U_inv << std::endl;  
  std::cout << "F2 = " << F2 << "     Fdiff = " << Fdiff << std::endl;
  std::cout << "Error of the self-consistent equation = " << scdiff << std::endl;  
    
  return true;
}

double negative_mean_field_function(vec const& x){
  double delta = x[0];
  bool verbose = false;
  double mu = sci2b.calc_chemical_potential(delta, verbose);
  sci2b.set_input(delta, mu);
  return - sci2b.calc_mean_field_function();
}

bool compare_mean_field_functions(vec const& x1, vec const& x2){
  double f1 = negative_mean_field_function(x1);
  double f2 = negative_mean_field_function(x2);
  return f1 < f2;
}

bool find_first_order_transition_point_bilayer2(rpa::parameters const& pr, double& U, double delta_i){
  std::cout << "Finding a value of U at the maximum value of the mean field function..." << std::endl;  

  /* Getting parameters */
  int L = pr.L;
  double T = pr.T;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Setting the parameters */
  double delta0 = 0.;
  double mu0 = 0.;
  bool non_zero_delta = true;  
  sci2b.set_parameters(pr, L, *ts, U, filling, T, delta0, mu0, continuous_k, non_zero_delta);

  /* Setting NelderMead class */
  int dim = 1;
  NelderMead nm(dim);
  nm.f = &negative_mean_field_function;
  nm.compare = &compare_mean_field_functions;
  nm.set_eps(1e-6);
  nm.reset();
  double delta_max = pr.find_U1st_delta_max;
  double delta_min = pr.find_U1st_delta_min;  
  if (delta_i > delta_max || delta_i < delta_min) {
    delta_i = 0.5 * (delta_max + delta_min);
  }
  vec x0{delta_i};
  vec x1{delta_max};
  std::vector<vec> xs(dim+1);
  xs[0] = x0;
  xs[1] = x1;
  nm.init_x(xs);
  nm.sort();

  /* Output */
  ofstream out_nm("optimal_delta_Nelder-Mead.text");
  out_nm << std::setprecision(12);
  
  /* Finding an optimal value. */
  std::size_t max_iter = 100;
  bool optimized = true;
  for(std::size_t t=0; t < max_iter; ++t){
    nm.output(out_nm);
    nm.step();
    if (nm.is_terminated()) { break; }      
    if (t == max_iter - 1){ optimized = false; }
  }
  nm.output(out_nm);

  /* Result */
  double f_opt = 0.;
  nm.get_result(x0, f_opt);
  double delta = x0[0];
  bool valid = delta >= delta_min && delta <= delta_max;  
  if (optimized && valid){    
    /* For the ordered state */
    bool verbose = false;
    double T = 0.0;
    double mu = sci2b.calc_chemical_potential(delta, verbose);
    sci2b.set_input(delta, mu);
    double F2 = sci2b.calc_energy();
    
    /* For the disordered state */
    double delta0 = 0;
    double mu0 = sci2b.calc_chemical_potential(delta0, verbose);
    sci2b.set_input(delta0, mu0);
    double F1 = sci2b.calc_energy();
    
    /* Comparing the energies. */
    out_nm << "F2 - F1 = " << F2 - F1 << std::endl;
    if (F2 < F1) {
      out_nm << "Optimal point: ( " << delta << " " << f_opt << " )" << std::endl;
      U = - 1. / f_opt;  // Negative      
      return true;
    } else {
      out_nm << "Not found." << std::endl;
      out_nm.close();    
      return false;
    }
  } else {
    std::cout << "Number of iteration reaches the limit " << max_iter << std::endl;
    std::cout << "The remained error squared is " << std::setprecision(12) << f_opt << std::endl;
    out_nm << "Not found." << std::endl;
    out_nm.close();    
    return false;
  }  
}

double calc_energy_diff_bilayer_anneal(rpa::parameters& pr, double U, double U_max, double U_delta){
  double F1 = 0., F2 = 0., delta = 0., mu = 0., ch_gap = 0., diff = 0.;
  bool set_init_val = false;

  /* Gradually decreasing U. */
  for(double Ua = U_max; Ua > U; Ua -= U_delta){
    std::cout << "Annealing Ua = " << Ua << "   to U = " << U << std::endl;      
    std::tie(F1, F2, delta, mu, ch_gap, diff) = calc_total_energies_bilayer(pr, Ua, delta, mu, set_init_val);
    set_init_val = true;
    std::cout << "delta = " << delta << "   mu = " << mu << "   charge_gap = " << ch_gap << "   sc_diff = " << diff << std::endl;
  }

  /* For U */
  std::tie(F1, F2, delta, mu, ch_gap, diff) = calc_total_energies_bilayer(pr, U, delta, mu, set_init_val);
  std::cout << "delta = " << delta << "   mu = " << mu << "   charge_gap = " << ch_gap << "   sc_diff = " << diff << std::endl;
  
  return F2 - F1;
}

double find_first_order_transition_point_bilayer_anneal(rpa::parameters& pr, double Uc){
  std::cout << "Finding a first-order transition point using an annealing process." << std::endl;

  /* Extracting parameters. */
  double U_max = pr.find_U1st_U_max;
  double U_min = pr.find_U1st_U_min;  
  double U_delta = pr.find_U1st_U_delta;  

  /* Target value */
  double target = 0.;

  /* Function */
  using std::placeholders::_1;
  auto eq = std::bind(calc_energy_diff_bilayer_anneal, std::ref(pr), _1, U_max, U_delta);

  BinarySearch bs(pr.continuous_k);
  bs.set_x_MIN(U_min);
  bs.set_x_MAX(U_max);

  /* Initial value */
  double U = U_max - U_delta;
  double Ui = Uc + 10. * U_delta;
  if (Ui < U) { U = Ui; }
  
  bool sol_found = bs.find_solution(U, target, eq);

  if ( sol_found ) {
    return U;
  } else {
    return 0;
  }  
}

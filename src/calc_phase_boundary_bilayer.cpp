/*****************************************************************************
*
* Functions for obtaining the phase boundary.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_phase_boundary.h"
#include "BinarySearch.h"
#include "find_critical_point.h"
#include "find_critical_U.h"
#include "find_metal_insulator_transition.h"
#include "self_consistent_eq.h"
#include "calc_chemical_potential.h"

/* Instances */
extern SelfConsistentIntegrand2Bilayer sci2b;

void calc_phase_boundary_U_bilayer(path& base_dir, rpa::parameters& pr){
  std::cout << "Obtaining the phase boundary as a function of U..." << std::endl;
    
  /* Output */
  ofstream out_pb;
  out_pb.open( base_dir / "phase-boundary-U.out");
  out_pb << "# U     tz_c" << std::endl;
  
  /* Precision */
  int prec = 12;
  
  /* U values */
  int n_U = int((pr.U_max - pr.U_min) / pr.U_delta + 1e-12) + 1;
  std::vector<double> Us(n_U);
  for(int Ui=0; Ui < n_U; Ui++){ Us[Ui] = pr.U_min + pr.U_delta * Ui; }

  /* For each U */
  for(int Ui=0; Ui < n_U; Ui++){
    double U = Us[Ui];
    std::cout << "U = " << U << std::endl;
    if ( std::abs(U) < 1e-12 ) {
      std::cout << "Skipping the case of U = 0." << std::endl;
      continue;
    }
    pr.U = U;
    if ( pr.fix_J ) {
      pr.t1 = 0.5*sqrt(pr.J*U);
    }

    if (pr.find_metal_insulator_transition) {
      double tz_c = find_metal_insulator_transition_t4_bilayer(pr);
      out_pb << std::setprecision(prec) << U << std::setw(prec+8) << tz_c << std::endl;      
    } else {
      double tz_c = find_critical_point_bilayer(pr);
      out_pb << std::setprecision(prec) << U << std::setw(prec+8) << tz_c << std::endl;
    }
  }
  
  out_pb.close();
}

double calc_energy_diff(rpa::parameters& pr, double U, bool anneal = false, double U_max = 0., double U_delta = 0.){
  double F1 = 0., F2 = 0., delta = 0., mu = 0., ch_gap = 0., diff = 0.;
  bool set_init_val = false;

  /* Gradually decreasing U. */
  if (anneal) {
    for(double Ua = U_max; Ua > U; Ua -= U_delta){
      std::cout << "Annealing Ua = " << Ua << "   to U = " << U << std::endl;      
      std::tie(F1, F2, delta, mu, ch_gap, diff) = calc_total_energies(pr, Ua, delta, mu, set_init_val);
      set_init_val = true;
    }
  }

  /* For U */
  std::tie(F1, F2, delta, mu, ch_gap, diff) = calc_total_energies(pr, U, delta, mu, set_init_val);
  
  return F2 - F1;
}

double calc_grand_potential_bilayer(SelfConsistentIntegrand2Bilayer& sc, double delta){
  double mu = sc.calc_chemical_potential(delta);
  sc.set_input(delta, mu);
  return sc.calc_energy();
}

double find_first_order_transition_point_bilayer(rpa::parameters& pr, double delta_i){
  std::cout << "Finding a first-order transition point." << std::endl;

  /* Getting parameters */
  int L = pr.L;
  double U = pr.U;
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
  double F0 = calc_grand_potential_bilayer(sci2b, delta0);
  std::cout << "F0 = " << F0 << std::endl;
  
  /* Finding the order parameter (delta) to match the disordered grand potential. */
  double target = F0;
  using std::placeholders::_1;
  auto eq = std::bind(calc_grand_potential_bilayer, std::ref(sci2b), _1);

  BinarySearch bs(continuous_k);
  bs.set_x_MIN(1e-5);
  bs.set_x_MAX(0.5);

  std::cout << "Binary search to find a solution." << std::endl;  
  double delta = delta_i;
  // bool sol_found = bs.find_solution(delta, target, eq, false, 0., true);
  bool sol_found = bs.find_solution(delta, target, eq);  
  if (!sol_found) {
    return 0.;
  }

  /* Calculating the corresponding value of U. */
  double mu = sci2b.calc_chemical_potential(delta);
  sci2b.set_input(delta, mu);
  double U_inv = sci2b.calc_mean_field_function();
  sci2b.set_U(1./U_inv);
  
  /* Output */
  double F1 = sci2b.calc_energy();
  double Fdiff = F1 - F0;
  double scdiff = sci2b.calc_diff();
  std::cout << "delta = " << delta << "   mu = " << mu << "   U_inv = " << U_inv << "   U1st = " << 1. / U_inv << std::endl;  
  std::cout << "F1 = " << F1 << "     Fdiff = " << Fdiff << std::endl;
  std::cout << "Error of the self-consistent equation = " << scdiff << std::endl;  
    
  return 1. / U_inv;
}

// double find_first_order_transition_point_old(rpa::parameters& pr, double Uc){
//   std::cout << "Finding a first-order transition point." << std::endl;

//   /* Extracting parameters. */
//   bool anneal = pr.find_U1st_anneal;
//   double U_max = pr.find_U1st_U_max;
//   double U_delta = pr.find_U1st_U_delta;  

//   /* Target value */
//   double target = 0.;

//   /* Function */
//   using std::placeholders::_1;
//   auto eq = std::bind(calc_energy_diff, std::ref(pr), _1, anneal, U_max, U_delta);

//   BinarySearch bs(pr.continuous_k);

//   /* The first-order transition point is expected to be greater than Uc. */
//   bs.set_x_MIN(Uc + 1e-12);
//   bs.set_x_MAX(U_max);

//   double U = Uc + 10. * U_delta;
//   bool sol_found = bs.find_solution(U, target, eq);

//   if ( sol_found ) {
//     return U;
//   } else {
//     return 0;
//   }  
// }

void calc_phase_boundary_t4_bilayer(path& base_dir, rpa::parameters& pr){
  std::cout << "Obtaining the phase boundary as a function of t4..." << std::endl;
    
  /* Output */
  ofstream out_pb;
  out_pb.open( base_dir / "phase-boundary-t4.out");
  out_pb << "# tz     U_c(possible)     U_1st" << std::endl;
  
  /* Precision */
  int prec = 12;
  int sw = prec + 10;
  
  /* t4 values */
  int n_t4 = int((pr.t4_max - pr.t4_min) / pr.t4_delta + 1e-12) + 1;
  std::vector<double> t4s(n_t4);
  for(int t4i=0; t4i < n_t4; t4i++){ t4s[t4i] = pr.t4_min + pr.t4_delta * t4i; }

  /* For each t4 */
  for(int t4i=0; t4i < n_t4; t4i++){
    double t4 = t4s[t4i];
    if ( std::abs(t4) < 1e-12 ) {
      std::cout << "Skipping the case of t4 = 0." << std::endl;
      continue;
    }

    /* Setting t4 */
    pr.t4 = t4;    
    if (pr.check_mean_field_function) {
      check_mean_field_eq_bilayer(base_dir, pr);
    }    

    /* Finding the possible critical point. */
    double Uc = find_critical_U_bilayer(pr);   // Using delta == 0    
    std::cout << "Uc = " << Uc << std::endl;
    
    /* Checking if a first-order transition occurs. */
    if (pr.find_first_order_transition) {
      if (pr.check_details) {
	/* Output */
	double shift = std::pow(10, 6);  
	double t4r = std::round(pr.t4*shift)/shift;	
	std::string ofn = "free_energy-t4_"+std::to_string(t4r)+".text";
	ofstream out_e(base_dir/ofn);
	out_e << "# tz     U   F_disorder   F_order   delta_order   mu_order   ch_gap   diff" << std::endl;
    
	std::cout << "Checking the energies for each U." << std::endl;
	double U_delta = 0.001;
	double U_max = Uc + U_delta * 60;
	double delta = 0., mu = 0., ch_gap = 0., diff = 0.;
	bool set_init_val = false;
	for(double U = U_max; U > Uc + 1e-12; U -= U_delta){
	  double F1 = 0., F2 = 0.;
	  std::tie(F1, F2, delta, mu, ch_gap, diff) = calc_total_energies(pr, U, delta, mu, set_init_val);
	  set_init_val = true;
	  out_e << std::setw(sw) << t4 << std::setw(sw) << U << std::setw(sw) << F1 << std::setw(sw) << F2 << std::setw(sw) << delta << std::setw(sw) << mu << std::setw(sw) << ch_gap << std::setw(sw) << diff << std::endl;
	}
	out_e.close();
      }

      double delta_i = pr.init_value;
      double U1st = find_first_order_transition_point_bilayer(pr, delta_i);
      
      /* Output */
      out_pb << std::setprecision(prec) << t4 << std::setw(sw) << Uc << std::setw(sw) << U1st << std::endl;
    } else {
      /* Output */
      out_pb << std::setprecision(prec) << t4 << std::setw(sw) << Uc << std::endl;
    }
  }

  out_pb.close();  
}

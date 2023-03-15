/*****************************************************************************
*
* Functions for finding a metal-insulator transition point.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "BinarySearch.h"
#include "calc_chemical_potential.h"
#include "find_metal_insulator_transition.h"

double calc_charge_gap_bilayer_t4(rpa::parameters& pr, double t4, double delta){
  /* Setting the parameter */
  pr.t4 = t4;

  /* Extracting parameters */
  int L = pr.L;

  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Calculating the charge gap */
  double ch_gap = 0., mu = 0.;
  bool negative = true;
  std::tie(ch_gap, mu) = calc_charge_gap_bilayer(L, *ts, delta, negative);

  return ch_gap;
}

double find_metal_insulator_transition_t4_bilayer(rpa::parameters& pr){
  std::cout << "Finding a metal-insulator transition point..." << std::endl;
  
  /* Asserts */
  assert(pr.filling == 0.5);  // Assume half filling.
  assert(pr.continuous_k == false);  // Assume a finite size system.

  /* Assume no order. */
  double delta = 0.;
  
  /* Setting the parameter to an initial value. */
  double tz = pr.init_value;
  
  /* Binary search */
  using std::placeholders::_1;
  auto scc = std::bind(calc_charge_gap_bilayer_t4, std::ref(pr), _1, delta);
  BinarySearch bs(pr.continuous_k);
  bs.set_x_MIN(0);
  double target = 0.;  
  bool verbose = true;
  bool additive = false;
  double x_delta = 0.;
  bool sol_found = bs.find_solution(tz, target, scc, additive, x_delta, verbose);
  
  double TP = 0;
  if ( sol_found ) {
    TP = tz;
  } else {
    std::cerr << "A solution was not found." << std::endl;
    TP = 0;
  }
  
  return TP;
}

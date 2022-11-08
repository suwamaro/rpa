/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "self_consistent_eq.h"
#include "rpa_util.h"
#include "calc_gap.h"
#include "calc_single_particle_energy.h"
#include "calc_chemical_potential.h"
#include "BinarySearch.h"
#include "find_critical_T.h"

double self_consistent_eq_two_site(hoppings_two_site const& ts, double T, double delta){  
  double k = 0.;
  cx_double ek1 = ts.ek1(0, 0, 0);  
  double zki_over_delta = zk_over_delta(ek1, ts.tz, k, delta);

  /* Fermi density */
  double Ep = eigenenergy_HF(1, ek1, 0, 0, ts.tz, k, delta);
  double Em = eigenenergy_HF(-1, ek1, 0, 0, ts.tz, k, delta);  
  double mu = 0.5 * ( Em + Ep );
  double n_minus = fermi_density(Em, kB*T, mu);
  double n_plus = fermi_density(Ep, kB*T, mu);	    

  /* The self-consistent condition: sum = 1/U */  
  return zki_over_delta * ( n_minus - n_plus ) / 2.0;
}

double solve_self_consistent_eq_two_site(hoppings_two_site const& ts, double U, double T){
  
  std::cout << "Finding a self-consistent solution for U = " << U << ", T = " << T << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_two_site, std::ref(ts), T, _1 );

  bool continuous_k = false;
  BinarySearch bs(continuous_k);
  bool sol_found = bs.find_solution( delta, target, scc );

  if ( sol_found ) {
    return delta;
  } else {
    return 0;
  }  
}

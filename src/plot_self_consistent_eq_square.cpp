/*****************************************************************************
*
* Functions for checking the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "plot_self_consistent_eq_square.h"
#include "self_consistent_eq.h"

void plot_self_consistent_eq_square(double U){
  /* Calculating the value of the self-consistent equation as a function of delta */
  std::cout << "U=" << U << std::endl;
  
  double t = 1.;
  double mu = 0;
  int L = 256;
  
  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_square, L, t, mu, _1 );

  std::ofstream out;
  out.open("sceq-delta.text");

  int prec = 15;
    
  double delta_delta = 0.001;
  double delta_max = 0.5 * U;
  for(double delta = delta_delta; delta <= delta_max; delta += delta_delta){
    double val = scc( delta );
    out << delta << std::setw( prec ) << val << std::endl;
  }
  double delta0 = solve_self_consistent_eq_square( L, t, mu, U );
  std::cout << "delta = " << delta0 << std::endl;

  out.close();
}

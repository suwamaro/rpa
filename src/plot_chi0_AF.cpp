/*****************************************************************************
*
* Functions for plotting the susceptibility with respect to the antiferromagnetic 
* ground state.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_gap.h"
#include "self_consistent_eq.h"

void plot_chi0_AF(double U){
  /* Calculating chi0 with respect to the AF ground state as a function of omega */
  double t = 1.;
  double mu = 0;
  int L = 48;
  double k1 = 2. * M_PI / L;
  
  double delta = solve_self_consistent_eq_square( L, t, mu, U );
  std::cout << "delta = " << delta << std::endl;
  
  boost::filesystem::ofstream out;
  out.open("chi0-omega.text");

  int prec = 15;
  double qx = M_PI;
  double qy = 0;
  
  double omega_delta = 0.0001;
  double omega_MIN = 0.1;
  double omega_MAX = 0.4;
  for(double omega = omega_MIN; omega <= omega_MAX; omega += omega_delta){
    double eigval = calc_eigval( L, t, mu, U, delta, qx, qy, omega );
    out << omega << std::setw( prec ) << eigval << std::endl;
  }
  double omega0 = calc_gap( L, t, mu, U, delta, qx, qy );
  
  std::cout << "omega0 = " << omega0 << std::endl;
  
  out.close();
}

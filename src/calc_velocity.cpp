/*****************************************************************************
*
* Functions for calculating the velocity
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_velocity.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"

void calc_velocity_square() {  
  /* Calculating the velocity */
  double t = 1.;
  double mu = 0;
  int L = 256;
  double k1 = 2. * M_PI / (double)L;
  double qx = M_PI;
  double qy = M_PI - k1;
  int prec = 15;
  
  boost::filesystem::ofstream out;
  out.open("velocity.text");
  
  double U_delta = 0.05;
  double U_min = 0.8;
  double U_max = 15.0;
  for(double U = U_min; U <= U_max; U += U_delta){
    double delta = solve_self_consistent_eq_square( L, t, mu, U );
    std::cout << "delta = " << delta << std::endl;
    double E1 = calc_gap_square( L, t, mu, U, delta, qx, qy );
    // double J = 4. * t * t / U;
    double velocity = E1 / k1;
    out << U << std::setw( prec ) << velocity << std::endl;
  }
  out.close();
}

void calc_velocity_cubic() {  
  /* Calculating the velocity */
  double t = 1.;
  double mu = 0;
  int L = 32;
  double k1 = 2. * M_PI / (double)L;
  double qx = M_PI;
  double qy = M_PI;
  double qz = M_PI - k1;
  int prec = 15;
  
  boost::filesystem::ofstream out;
  out.open("velocity.text");
  
  double U_delta = 1.;
  double U_min = 1.;
  double U_max = 2000.0;
  for(double U = U_min; U <= U_max; U += U_delta){
    double delta = solve_self_consistent_eq_cubic( L, t, mu, U );
    std::cout << "delta = " << delta << std::endl;
    double E1 = calc_gap_cubic( L, t, mu, U, delta, qx, qy, qz );
    // double J = 4. * t * t / U;
    double velocity = E1 / k1;
    out << U << std::setw( prec ) << velocity << std::endl;
  }
  out.close();
}

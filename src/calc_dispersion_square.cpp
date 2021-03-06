/*****************************************************************************
*
* Functions for calculating the dispersion
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_dispersion.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"
#include "BinarySearch.h"

void calc_dispersion_square(double U, int L){
  double t = 1.;
  double t_bar = 0;
  double mu = 0;
  double k1 = 2. * M_PI / (double)L;

  std::unique_ptr<hoppings_square> ts;
  ts = hoppings_square::mk_square(t, t_bar);
  
  double delta = solve_self_consistent_eq_square( L, *ts, mu, U );
  std::cout << "delta = " << delta << std::endl;
  
  std::ofstream out;
  out.open("dispersion.text");

  int prec = 15;

  double qx = 0;
  double qy = 0;
  int q_idx = 0;
  
  /* Square lattice */
  /* M -> X -> Gamma -> M -> R -> X */
  /* Gamma = ( 0, 0 )               */
  /* M = ( pi, 0 )                  */
  /* X = ( pi/2, pi/2 )             */
  /* R = ( pi, pi )                 */

  auto output_omega = [&](){
    std::cout << "( qx, qy ) = ( " << qx << ", " << qy << " )" << std::endl;

    /* Gap is 0 at (0,0) and (pi,pi) */
    if ( std::abs( qx ) + std::abs( qy ) < 1e-12 ) {
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << 0 << std::setw( prec ) << U << std::endl;
    } else if ( std::abs( qx - M_PI ) + std::abs( qy - M_PI ) < 1e-12 ) {
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << 0 << std::setw( prec ) << U << std::endl;
    } else {
      /* Finding the pole of the RPA susceptibility. */
      double omega = calc_gap_square( L, t, mu, U, delta, qx, qy );
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omega << std::setw( prec ) << U << std::endl;
    }
  };
      
  for(int x=0; x < L/4; x++){
    qx = M_PI - k1 * x;
    qy = k1 * x;
    output_omega();
    ++q_idx;
  }

  for(int x=0; x < L/4; x++){
    qx = 0.5 * M_PI - k1 * x;
    qy = 0.5 * M_PI - k1 * x;
    output_omega();
    ++q_idx;
  }
  
  for(int x=0; x < L/2; x++){
    qx = k1 * x;
    qy = 0;
    output_omega();
    ++q_idx;
  }

  for(int y=0; y < L/2; y++){
    qx = M_PI;
    qy = k1 * y;
    output_omega();
    ++q_idx;
  }

  for(int x=0; x <= L/4; x++){
    qx = M_PI - k1 * x;
    qy = M_PI - k1 * x;
    output_omega();
    ++q_idx;
  }

  out.close();
}

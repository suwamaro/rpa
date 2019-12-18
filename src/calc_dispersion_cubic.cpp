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

void calc_dispersion_cubic(double U, int L){
  double t = 1.;
  double mu = 0;
  double k1 = 2. * M_PI / (double)L;

  double delta = solve_self_consistent_eq_cubic( L, t, mu, U );
  std::cout << "delta = " << delta << std::endl;
  
  std::ofstream out;
  out.open("dispersion.text");

  int prec = 15;

  double qx = 0;
  double qy = 0;
  double qz = 0;
  int q_idx = 0;
  
  /* Simple cubic lattice */
  /* Gamma -> X -> M -> R -> Gamma */
  /* Gamma = ( 0, 0, 0 )           */
  /* X = ( pi, 0, 0 )              */
  /* M = ( pi, pi, 0 )             */
  /* R = ( pi, pi, pi )            */

  auto output_omega = [&](){
    std::cout << "( qx, qy, qz ) = ( " << qx << ", " << qy << ", " << qz << " )" << std::endl;

    /* Gap is 0 at (0,0,0) and (pi,pi,pi) */
    if ( std::abs( qx ) + std::abs( qy ) + std::abs( qz ) < 1e-12 ) {
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << 0 << std::setw( prec ) << U << std::endl;
    } else if ( std::abs( qx - M_PI ) + std::abs( qy - M_PI ) + std::abs( qz - M_PI ) < 1e-12 ) {
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << 0 << std::setw( prec ) << U << std::endl;
    } else {
      /* Finding the pole of the RPA susceptibility. */
      double omega = calc_gap_cubic( L, t, mu, U, delta, qx, qy, qz );
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << omega << std::setw( prec ) << U << std::endl;
    }
  };
      
  for(int x=0; x < L/2; x++){
    qx = k1 * x;
    qy = 0;
    qz = 0;
    output_omega();
    ++q_idx;
  }

  for(int y=0; y < L/2; y++){
    qx = M_PI;
    qy = k1 * y;
    qz = 0;
    output_omega();
    ++q_idx;
  }
  
  for(int z=0; z < L/2; z++){
    qx = M_PI;
    qy = M_PI;
    qz = k1 * z;
    output_omega();
    ++q_idx;
  }

  for(int x=0; x <= L/2; x++){
    qx = M_PI - k1 * x;
    qy = qx;
    qz = qx;
    output_omega();
    ++q_idx;
  }

  out.close();
}

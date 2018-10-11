/*****************************************************************************
*
* Functions for calculating the spectrum
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_spectrum.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"

void calc_spectrum_square(double U){
  /* Parameters */
  double t = 1.;
  double mu = 0;
  int L = 16;
  double eta = 0.001;
  double k1 = 2. * M_PI / (double)L;
  // double delta_omega = 0.05;
  // double max_omega = 0.001;

  double delta_omega = 0.01;
  double max_omega = 2.0;
  
  int prec = 15;
  
  /* Omegas */
  int n_omegas = int(max_omega/delta_omega+1);
  std::vector<double> omegas(n_omegas);
  for(int o=0; o < n_omegas; o++){ omegas[o] = delta_omega * o; }

  /* Calculate the gap */
  double delta = solve_self_consistent_eq_square( L, t, mu, U );
  std::cout << "delta = " << delta << std::endl;

  /* Output */
  boost::filesystem::ofstream out;
  out.open("spectrum.text");

  /* Wavenumbers */
  double qx = 0;
  double qy = 0;
  int q_idx = 0;
  
  /* Square lattice */
  /* M -> X -> Gamma -> M -> R -> X */
  /* Gamma = ( 0, 0 )               */
  /* M = ( pi, 0 )                  */
  /* X = ( pi/2, pi/2 )             */
  /* R = ( pi, pi )                 */

  auto output_spectrum = [&](){
    std::cout << "( qx, qy ) = ( " << qx << ", " << qy << " )" << std::endl;

    /* Gap is 0 at (0,0) and (pi,pi) */
    if ( std::abs( qx ) + std::abs( qy ) < 1e-12 ) {
      for(int o=0; o < n_omegas; o++){
	out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omegas[o] << std::setw( prec ) << 0 << std::setw( prec ) << U << std::endl;
      }
    } else if ( std::abs( qx - M_PI ) + std::abs( qy - M_PI ) < 1e-12 ) {
      for(int o=0; o < n_omegas; o++){
	out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omegas[o] << std::setw( prec ) << 0 << std::setw( prec ) << U << std::endl;
      }
    } else {
      /* Finding the pole of the RPA susceptibility. */
      for(int o=0; o < n_omegas; o++){
	cx_double cx_omega(omegas[o], eta);
	double chi = calc_intensity_square( L, t, mu, U, delta, qx, qy, cx_omega );
	out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omegas[o] << std::setw( prec ) << chi << std::setw( prec ) << U << std::endl;
      }
    }
  };
      
  for(int x=0; x < L/4; x++){
    qx = M_PI - k1 * x;
    qy = k1 * x;
    output_spectrum();
    ++q_idx;
  }

  for(int x=0; x < L/4; x++){
    qx = 0.5 * M_PI - k1 * x;
    qy = 0.5 * M_PI - k1 * x;
    output_spectrum();
    ++q_idx;
  }
  
  for(int x=0; x < L/2; x++){
    qx = k1 * x;
    qy = 0;
    output_spectrum();
    ++q_idx;
  }

  for(int y=0; y < L/2; y++){
    qx = M_PI;
    qy = k1 * y;
    output_spectrum();
    ++q_idx;
  }

  for(int x=0; x <= L/4; x++){
    qx = M_PI - k1 * x;
    qy = M_PI - k1 * x;
    output_spectrum();
    ++q_idx;
  }

  out.close();
}

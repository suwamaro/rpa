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

void calc_spectrum_square(double U, int L, double eta){
  /* Parameters */
  double t = 1.;
  double mu = 0;
  double k1 = 2. * M_PI / (double)L;
  
  int prec = 15;
  
  /* Omegas */
  double delta_omega = 0.01;
  // double delta_omega = 0.002;
  
  double max_omega = 10.0;
  int n_omegas = int(max_omega/delta_omega+0.5);
  std::vector<double> omegas(n_omegas);
  for(int o=1; o <= n_omegas; o++){ omegas[o-1] = delta_omega * o; }

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
    std::cout << "( qidx, qx, qy ) = ( " << q_idx << ", " << qx << ", " << qy << " )" << std::endl;
    for(int o=0; o < n_omegas; o++){
      cx_double cx_omega(omegas[o], eta);

      /* Taking the linear combitation from 0 (index) at (qx,qy) and from 1 at (qx+M_PI,qy+M_PI) */
      double chi = 0;
      chi += calc_intensity_square( L, t, mu, U, delta, qx, qy, cx_omega, 0 );
      chi += calc_intensity_square( L, t, mu, U, delta, qx + M_PI, qy + M_PI, cx_omega, 1 );
      chi *= 0.5;

      /* Output */
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omegas[o] << std::setw( prec ) << chi << std::setw( prec ) << U << std::endl;
    }
  };

  // // for check
  // qx = k1;
  // qy = 0;
  // output_spectrum();
  // return;
  
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

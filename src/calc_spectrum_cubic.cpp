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

void calc_spectrum_cubic(double U, int L, double eta){
  /* Parameters */
  double t = 1.;
  double mu = 0;
  double k1 = 2. * M_PI / (double)L;
  
  int prec = 15;
  
  /* Omegas */
  // double delta_omega = 0.01;
  double delta_omega = 0.001;
  
  // double max_omega = 4.0;
  double max_omega = 10.0;
  int n_omegas = int(max_omega/delta_omega+0.5);
  std::vector<double> omegas(n_omegas);
  for(int o=1; o <= n_omegas; o++){ omegas[o-1] = delta_omega * o; }

  /* Calculate the gap */
  double delta = solve_self_consistent_eq_cubic( L, t, mu, U );
  std::cout << "delta = " << delta << std::endl;

  /* Output */
  boost::filesystem::ofstream out;
  out.open("spectrum.text");

  /* Wavenumbers */
  double qx = 0;
  double qy = 0;
  double qz = 0;
  int q_idx = 0;
  
  /* Cubic lattice */
  /* M -> X -> Gamma -> M -> R -> X */
  /* Gamma = ( 0, 0 )               */
  /* M = ( pi, 0 )                  */
  /* X = ( pi/2, pi/2 )             */
  /* R = ( pi, pi )                 */

  /* From the response function to the dynamic structure factor */
  double factor_dsf = 2.;
  
  auto output_spectrum = [&](){
    std::cout << "( qidx, qx, qy, qz ) = ( " << q_idx << ", " << qx << ", " << qy << ", " << qz << " )" << std::endl;
    for(int o=0; o < n_omegas; o++){
      cx_double cx_omega(omegas[o], eta);

      /* The (S^+ S^-) response function, or the retarded Green's function */
      cx_double chi = calc_intensity_cubic( L, t, mu, U, delta, qx, qy, qz, cx_omega);
      double spec = factor_dsf * std::imag(chi);
      
      /* Output */
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << omegas[o] << std::setw( prec ) << spec << std::setw( prec ) << U << std::endl;
    }
  };

  // // for check with the bilayer
  // for(int z=0; z < 2; z++){
  //   qz = M_PI * z;
    
  //   for(int x=0; x < L/4; x++){
  //     qx = M_PI - k1 * x;
  //     qy = k1 * x;
  //     output_spectrum();
  //     ++q_idx;
  //   }
    
  //   for(int x=0; x < L/4; x++){
  //     qx = 0.5 * M_PI - k1 * x;
  //     qy = 0.5 * M_PI - k1 * x;
  //     output_spectrum();
  //     ++q_idx;
  //   }
    
  //   for(int x=0; x < L/2; x++){
  //     qx = k1 * x;
  //     qy = 0;
  //     output_spectrum();
  //     ++q_idx;
  //   }
    
  //   for(int y=0; y < L/2; y++){
  //     qx = M_PI;
  //     qy = k1 * y;
  //     output_spectrum();
  //     ++q_idx;
  //   }
    
  //   for(int x=0; x <= L/4; x++){
  //     qx = M_PI - k1 * x;
  //     qy = M_PI - k1 * x;
  //     output_spectrum();
  //     ++q_idx;
  //   }
  // }
  
  for(int z=-1; z < 1; z++){
    qz = M_PI * z;
  
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
  }

  out.close();
}

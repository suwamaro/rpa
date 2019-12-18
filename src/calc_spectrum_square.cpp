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
  hoppings_simple ts(t);
  double k1 = 2. * M_PI / (double)L;
  
  int prec = 15;
  
  /* Omegas */
  double delta_omega = 0.001;
  double max_omega = ts.t_max() * 10;
  int n_omegas = int(max_omega/delta_omega+0.5);
  std::vector<double> omegas(n_omegas);
  for(int o=1; o <= n_omegas; o++){ omegas[o-1] = delta_omega * o; }

  /* Calculate the gap */
  double delta = solve_self_consistent_eq_square( L, t, mu, U );
  std::cout << "delta = " << delta << std::endl;

  /* Output */
  std::ofstream out_xy, out_z;
  out_xy.open("spectrum-xy.text");
  out_z.open("spectrum-z.text");  

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

  /* From the response function to the dynamic structure factor */
  double factor_dsf = 2.;
  
  auto output_spectrum = [&](){
    std::cout << "( qidx, qx, qy ) = ( " << q_idx << ", " << qx << ", " << qy << " )" << std::endl;
    for(int o=0; o < n_omegas; o++){
      cx_double cx_omega(omegas[o], eta);

      /* The (S^+ S^-) response function, or the retarded Green's function */
      double factor_dsf = 2.;
      cx_double chi_xy = calc_intensity_square( L, ts, mu, U, delta, qx, qy, cx_omega, false );
      double spec_xy = factor_dsf * std::imag(chi_xy);
      
      /* Output */
      out_xy << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omegas[o] << std::setw( prec ) << spec_xy << std::setw( prec ) << U << std::endl;

      /* The (S^z S^z) response function, or the retarded Green's function */
      cx_double chi_z = calc_intensity_square( L, ts, mu, U, delta, qx, qy, cx_omega, true );
      double spec_z = factor_dsf * std::imag(chi_z);
      
      /* Output */
      out_z << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omegas[o] << std::setw( prec ) << spec_z << std::setw( prec ) << U << std::endl;      
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

  out_xy.close();
  out_z.close();
}

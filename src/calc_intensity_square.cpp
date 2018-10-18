/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_gap.h"
#include "rpa_util.h"
#include <armadillo>

double calc_intensity_square(int L, double t, double mu, double U, double delta, double qx, double qy, cx_double omega){
  double k1 = 2. * M_PI / (double)L;
  
  cx_double A = 0, B = 0, D = 0;
  /* Summing up at all the wavevectors */
  for(int x=-L/2; x < L/2; x++){
    double kx = k1 * x;
    for(int y=-L/2; y < L/2; y++){
      double ky = k1 * y;
      
      double e_free = energy_free_electron( t, mu, kx, ky );
      double e_eps = 1e-12;
      if ( e_free < e_eps ) {
	/* Summing up over all k inside the Brillouin zone. */
	double diff_x = wave_vector_in_BZ( kx - qx );
	double diff_y = wave_vector_in_BZ( ky - qy );	
	double e_free2 = energy_free_electron( t, mu, diff_x, diff_y );
	add_to_sus_mat( A, B, D, e_free, e_free2, delta, omega );
      }
    }
  }

  int n_sites = L * L;
  A *= 2. / (double)n_sites;
  B *= 2. / (double)n_sites;
  D *= 2. / (double)n_sites;
  
  /* Prefactors */
  double factor_SS = 3.;

  /* From the response function to the dynamic structure factor */
  double factor_dsf = 2.;

  /* RPA */
  arma::cx_mat chi0_mat(2,2);
  chi0_mat(0,0) = A;
  chi0_mat(0,1) = B;
  chi0_mat(1,0) = std::conj(B); // Taking the conjugate  
  chi0_mat(1,1) = D;
  arma::cx_mat denom = arma::eye<arma::cx_mat>(2,2) - U * chi0_mat;
  arma::cx_mat chi_mat = chi0_mat * arma::inv(denom);

  /* Taking the linear combitation */
  cx_double chi = chi_mat(0,0) + chi_mat(0,1) + chi_mat(1,0) + chi_mat(1,1);

  /* Cancelling double counts related to the sublattices */
  chi *= 0.5;

  //   // for check
  // std::cout << qx << "   " << qy << "   " << omega << "   " << A << "   " << B << "   " << D << "   " << "   " << chi << std::endl;
  
  double chi_Im = std::imag( chi );
  if ( chi_Im > 0 ) {
    return factor_dsf * factor_SS * chi_Im;
  } else {   // Numerical error?
    return 0;
  }
  
}

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
  
  arma::cx_mat chi0(2,2);
  chi0(0,0) = A;
  chi0(0,1) = B;
  chi0(1,0) = std::conj(B);
  chi0(1,1) = D;
  arma::cx_mat denom = arma::eye<arma::cx_mat>(2,2) - U * chi0;
  arma::cx_mat chi = chi0 * arma::inv(denom);
  arma::vec eigval;
  arma::eig_sym( eigval, chi );
  // arma::cx_mat eigvec;
  // arma::eig_sym( eigval, eigvec, chi );
  
  // std::cout << eigval(0) << std::endl;
  // std::cout << eigval(1) << std::endl;
  // return arma::max(arma::abs(eigval));
  return arma::min(eigval);
  // return arma::max(eigval);
  
  
}

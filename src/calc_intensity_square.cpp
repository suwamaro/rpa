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

cx_double calc_intensity_square(int L, double t, double mu, double U, double delta, double qx, double qy, cx_double omega){
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
  
  /* RPA */
  arma::cx_mat chi0_mat(2,2);
  chi0_mat(0,0) = A;   // (A, A) correlation
  chi0_mat(0,1) = B;   // (A, B)
  chi0_mat(1,0) = B;   // (B, A)
  chi0_mat(1,1) = D;   // (B, B)
  arma::cx_mat denom = arma::eye<arma::cx_mat>(2,2) - U * chi0_mat;
  arma::cx_mat chi_mat = chi0_mat * arma::inv(denom);

  // Double counting from A and B
  double factor_sublattice = 0.5;
  cx_double chi = factor_sublattice * factor_sublattice * ( chi_mat(0,0) - chi_mat(1,0) - chi_mat(0,1) + chi_mat(1,1) );

  // // for check
  // std::cerr << "chi0: " << chi0_mat(0,0) << "  " << chi0_mat(0,1) << "  " << chi0_mat(1,0) << "  " << chi0_mat(1,1) << std::endl;

  // std::cerr << "chi:  " << chi_mat(0,0) << "  " << chi_mat(0,1) << "  " << chi_mat(1,0) << "  " << chi_mat(1,1) << std::endl;  

  // std::cerr << "return: " << chi << std::endl;
  
  return chi;  
}

// cx_double calc_intensity_square2(int L, double t, double mu, double U, double delta, double qx, double qy, cx_double omega, bool zz){  
//   double k1 = 2. * M_PI / (double)L;
//   double e_eps = 1e-12;  
//   cx_double polarization = 0;
  
//   /* Summing up over all the wavevectors in the magnetic BZ. */
//   for(int x=-L/2; x < L/2; x++){    
//     double kx = k1 * x;
    
//     for(int y=-L/2; y < L/2; y++){
//       double ky = k1 * y;
      
//       double e_free = energy_free_electron( t, mu, kx, ky );
//       if ( e_free < mu + e_eps ) {
// 	double factor = 1.;
// 	if ( std::abs( e_free - mu ) < e_eps ) {
// 	  factor = 0.5;
// 	}
	
// 	double kx2 = wave_vector_in_BZ( kx + qx );
// 	double ky2 = wave_vector_in_BZ( ky + qy );
// 	double e_free2 = energy_free_electron( t, mu, kx2, ky2 );

// 	double Ek = eigenenergy_HF_plus(e_free, delta);
// 	double Ek2 = eigenenergy_HF_plus(e_free2, delta);
	
// 	double element = 0;
// 	if ( zz ) {
// 	  element = 1. - (e_free * e_free2 + delta * delta) / (Ek * Ek2);
// 	} else {
// 	  element = 1. - (e_free * e_free2 - delta * delta) / (Ek * Ek2);
// 	}
	
// 	cx_double denom = 1. / (omega - Ek2 - Ek) + 1. / ( - std::conj(omega) - Ek2 - Ek);
// 	polarization += - factor * element / denom;
//       }
//     } /* end for y */
//   } /* end for x */

//   int nsites = L * L;
//   polarization /= (double)(2*nsites);
//   return polarization / ( 1. - U * polarization );
// }

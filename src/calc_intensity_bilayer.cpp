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

cx_double calc_intensity_bilayer(int L, hoppings const& ts, double mu, double U, double delta, double qx, double qy, double qz, cx_double omega, bool zz){
  double k1 = 2. * M_PI / (double)L;
  
  cx_double A = 0, B = 0, C = 0, D = 0;
  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){    
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;
	add_to_sus_mat2( ts, mu, A, B, C, D, qx, qy, qz, kx, ky, kz, delta, omega, zz );	
      }
    }
  }

  int n_sites = L * L * 2;
  double norm = 2. / (double)n_sites;
  A *= norm;
  B *= norm;
  C *= norm;
  D *= norm;
  
  /* RPA */
  arma::cx_mat chi0_mat(2,2);
  chi0_mat(0,0) = A;   // (A, A) correlation
  chi0_mat(0,1) = B;   // (A, B)
  chi0_mat(1,0) = C;   // (B, A)
  chi0_mat(1,1) = D;   // (B, B)

  /* Transverse = < \sigma^- \sigma^+ >; Longitudinal (zz) = < \sigma^z \sigma^z > */
  /* Note that 2 < \sigma^- \sigma^+ > = < \sigma^z \sigma^z > (U -> 0 for the SU(2) case) */  
  double factor_channel = 1.0;
  if ( zz ) { factor_channel = 0.5; }
  
  arma::cx_mat denom = arma::eye<arma::cx_mat>(2,2) - factor_channel * U * chi0_mat;
  arma::cx_mat chi_mat = chi0_mat * arma::inv(denom);
  
  // sigma-to-spin factor
  double factor_operator = 0.5;  
  cx_double chi = factor_operator * factor_operator * ( chi_mat(0,0) - chi_mat(1,0) - chi_mat(0,1) + chi_mat(1,1) );  
  
  return chi;  
}

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

cx_double calc_intensity_bilayer(int L, hoppings const& ts, double mu, double U, double delta, double qx, double qy, double qz, cx_double omega, int index){
  double k1 = 2. * M_PI / (double)L;
  
  cx_double A = 0, B = 0, C = 0, D = 0;

  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){    
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;
	add_to_sus_mat2( ts, mu, A, B, C, D, qx, qy, qz, kx, ky, kz, delta, omega );
	// add_to_sus_mat2( ts, mu, A, B, D, qx, qy, qz, kx, ky, kz, delta, omega );
	
      }
    }
  }

  int n_sites = L * L * 2;  
  A *= 2. / (double)n_sites;
  B *= 2. / (double)n_sites;
  C *= 2. / (double)n_sites;  
  D *= 2. / (double)n_sites;  

  /* RPA */
  arma::cx_mat chi0_mat(2,2);
  chi0_mat(0,0) = A;   // (A, A) correlation
  chi0_mat(0,1) = B;   // (A, B)
  chi0_mat(1,0) = C;   // (B, A)
  // chi0_mat(1,0) = B;   // (B, A)  
  // chi0_mat(1,0) = std::conj(B); // Taking the conjugate  
  chi0_mat(1,1) = D;   // (B, B)
  arma::cx_mat denom = arma::eye<arma::cx_mat>(2,2) - U * chi0_mat;
  arma::cx_mat chi_mat = chi0_mat * arma::inv(denom);

  // Double counting from A and B
  // Double counting in summing up for wavevectors because of the sublattice order?
  double factor_sublattice = 0.5;  
  // cx_double chi = factor_sublattice * factor_sublattice * ( chi_mat(0,0) - std::conj(chi_mat(1,0)) - std::conj(chi_mat(0,1)) + chi_mat(1,1) );
  // cx_double chi = factor_sublattice * factor_sublattice * ( chi_mat(0,0) - chi_mat(1,0) - std::conj(chi_mat(0,1)) + chi_mat(1,1) );
  
  cx_double chi = factor_sublattice * factor_sublattice * ( chi_mat(0,0) - chi_mat(1,0) - chi_mat(0,1) + chi_mat(1,1) );
  
  
  
  // // for check
  // std::cout << qx << "  " << qy << "  " << omega << "  " << A << "  " << B << "  " << D << "  " << chi_mat(0,0) << "  " << chi_mat(0,1) << "  " << chi_mat(1,0) << "  " << chi_mat(1,1) << "  " << chi << std::endl;

  return chi;  
}

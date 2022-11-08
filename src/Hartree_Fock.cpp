/*****************************************************************************
*
* Functions for the Hartree-Fock approximation.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "Hartree_Fock.h"

/* Band index */
int valence_index = 1;
int conduction_index = 0;

double wave_vector_in_BZ(double k){
  /* Returning - M_PI <= k < M_PI */
  double k_eps = 1e-12;
  while( k > M_PI - k_eps ) k -= 2. * M_PI;
  while( k < - M_PI - k_eps ) k += 2. * M_PI;
  if ( std::abs( k ) < k_eps ) k = 0;
  if ( std::abs( k + M_PI ) < k_eps ) k = - M_PI;
  return k;
}

double BZ_factor(double kx, double ky){  
  double mu_free = 0;  /* Assume at half filling */
  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc: t=1 */
  double factor = 0.0;
  double eps = 1e-12;  
  if ( e_free > mu_free + eps ) { factor = 0.0; }
  else if ( std::abs(e_free - mu_free) < eps ) { factor = 0.5; } /* On the zone boundary */
  else { factor = 1.0; }
  return factor;
}

int sublattice_spin_index(int g, int sigma){
  int idx = g << 1;
  if ( sigma == down_spin ) {
    idx += 1;
  }
  return idx;
}

int sign_spin_index(int s, int sigma){
  int idx = - s + 1;
  if ( sigma == down_spin ) {
    idx += 1;
  }
  return idx;
}

int band_spin_index(int b, int sigma){
  return (b << 1) | spin_index(sigma);
}

int spin_index(int sigma){
  return sigma == down_spin ? 1 : 0;
}

int index_to_spin(int idx){
  return idx == 0 ? up_spin : down_spin;
}

int index_to_sign(int idx){
  return idx == 0 ? 1 : -1;
}

cx_vec gs_HF1(int spin, int sign, cx_double ek1, cx_double tz, double kz, double delta){
  /* In the case where the Hamiltonian is block-diagonalized for each spin in the presence of the U(1) symmetry. */
  cx_double xki = xk(spin, ek1, tz, kz, delta);
  double zki = zk(ek1, tz, kz, delta);
  cx_double coef1 = xki*sqrt(0.5*(1+sign*spin*zki));
  cx_double coef2 = sign*sqrt(0.5*(1-sign*spin*zki));  
  
  cx_vec gs(2*NSUBL, arma::fill::zeros);
  if ( spin == up_spin ) {
    gs(0) = coef1;
    gs(2) = coef2;    
  } else {
    gs(1) = coef1;
    gs(3) = coef2;
  }
  return gs;
}

cx_mat gs_HF(cx_double ek1, cx_double tz, double kz, double delta){
  cx_vec HF_up_plus = gs_HF1(up_spin, 1, ek1, tz, kz, delta);
  cx_vec HF_up_minus = gs_HF1(up_spin, -1, ek1, tz, kz, delta);
  cx_vec HF_down_plus = gs_HF1(down_spin, 1, ek1, tz, kz, delta);
  cx_vec HF_down_minus = gs_HF1(down_spin, -1, ek1, tz, kz, delta);    
  
  cx_mat U(2*NSUBL, 2*NSUBL);

  /* Copying each vector to each column. */
  std::copy(HF_up_plus.begin(), HF_up_plus.end(), U.colptr(0));
  std::copy(HF_down_plus.begin(), HF_down_plus.end(), U.colptr(1));
  std::copy(HF_up_minus.begin(), HF_up_minus.end(), U.colptr(2));
  std::copy(HF_down_minus.begin(), HF_down_minus.end(), U.colptr(3));  
  
  return U;
}

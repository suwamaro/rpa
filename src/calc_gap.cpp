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

/* Reference: A. Singh and Z. Tesanovic, PRB 41, 11457 (1990) */

bool diagonal_H(double ek1, double ek2) {  
  return ek1 * ek1 + ek2 * ek2 < 1e-12;
}
double diagonal(double ek3, double delta, double E){
  return E - ek3 - delta;
}
cx_double off_diagonal(double ek1, double ek2){
  cx_double comp(ek1, - ek2);
  return comp;
}
cx_double calc_bk_up_in_minus(double ek1, double ek2, double ek3, double delta){
  double E = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double diag = diagonal(ek3, delta, E);
  cx_double off_diag = off_diagonal(ek1, ek2);
  double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
  return off_diag / denom;
  // return - off_diag / denom;  
}
cx_double calc_ak_up_in_minus(double ek1, double ek2, double ek3, double delta){
  // cx_double bk = calc_bk_up_in_minus(ek1, ek2, ek3, delta);
  // return sqrt(1. - std::norm(bk));  
  double E = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double diag = diagonal(ek3, delta, E);
  cx_double off_diag = off_diagonal(ek1, ek2);
  double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
  // return std::abs(diag) / denom;   /* Taking the absolute. */
  // return - diag / denom;
  return diag / denom;  
}
cx_double calc_ak_down_in_minus(double ek1, double ek2, double ek3, double delta){
  return calc_bk_up_in_minus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_down_in_minus(double ek1, double ek2, double ek3, double delta){
  return calc_ak_up_in_minus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_up_in_plus(double ek1, double ek2, double ek3, double delta){
  if ( diagonal_H(ek1, ek2) ) {
    return 1.;  /* It doesn't matter. */
  } else {  
    double E = eigenenergy_HF_plus(ek1, ek2, ek3, delta);
    double diag = diagonal(ek3, delta, E);
    cx_double off_diag = off_diagonal(ek1, ek2);
    double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
    return - off_diag / denom;
    // return - std::conj(off_diag) / denom;    /* for the square lattice */
    // return std::conj(off_diag) / denom;        
    // return off_diag / denom;    
  }
}
cx_double calc_ak_up_in_plus(double ek1, double ek2, double ek3, double delta){
  // cx_double bk = calc_bk_up_in_plus(ek1, ek2, ek3, delta);
  // return sqrt(1. - std::norm(bk));  
  if ( diagonal_H(ek1, ek2) ) {
    return 0.;
  } else {  
    double E = eigenenergy_HF_plus(ek1, ek2, ek3, delta);
    double diag = diagonal(ek3, delta, E);
    cx_double off_diag = off_diagonal(ek1, ek2);
    double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
    // return std::abs(diag) / denom;   /* Taking the absolute. */
    // return - diag / denom;
    return diag / denom;    
  }
}
cx_double calc_ak_down_in_plus(double ek1, double ek2, double ek3, double delta){
  return calc_bk_up_in_plus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_down_in_plus(double ek1, double ek2, double ek3, double delta){
  return calc_ak_up_in_plus(ek1, ek2, ek3, delta);
}

// double G_AA(double E, double delta, double spin){
//   return 0.5 * ( 1. - spin * delta / E );
// }
// double G_AA_minus(double ek1, double ek2, double ek3, double delta, double spin){
//   double E = ek3 - sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_AA( E, delta, spin );
// }
// double G_AA_plus(double ek1, double ek2, double ek3, double delta, double spin){
//   double E = ek3 + sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_AA( E, delta, spin );
// }
// double G_AA_up_minus(double ek1, double ek2, double ek3, double delta){
//   return G_AA_minus( ek1, ek2, ek3, delta, 1.0 );
// }
// double G_AA_down_minus(double ek1, double ek2, double ek3, double delta){
//   return G_AA_minus( ek1, ek2, ek3, delta, - 1.0 );
// }
// double G_AA_up_plus(double ek1, double ek2, double ek3, double delta){
//   return G_AA_plus( ek1, ek2, ek3, delta, 1.0 );
// }
// double G_AA_down_plus(double ek1, double ek2, double ek3, double delta){
//   return G_AA_plus( ek1, ek2, ek3, delta, - 1.0 );
// }

// cx_double G_AB(cx_double phi, double E){
//   return - 0.5 * phi / E;
// }
// cx_double G_AB_up_minus(double ek1, double ek2, double ek3, double delta){
//   cx_double phi(ek1, - ek2);  
//   double E = ek3 - sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_AB( phi, E );
// }
// cx_double G_AB_down_minus(double ek1, double ek2, double ek3, double delta){
//   cx_double phi(ek1, ek2);  
//   double E = ek3 - sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_AB( phi, E );
// }
// cx_double G_AB_up_plus(double ek1, double ek2, double ek3, double delta){
//   cx_double phi(ek1, - ek2);  
//   double E = ek3 + sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_AB( phi, E );
// }
// cx_double G_AB_down_plus(double ek1, double ek2, double ek3, double delta){
//   cx_double phi(ek1, ek2);  
//   double E = ek3 + sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_AB( phi, E );
// }

// cx_double G_BA_up_minus(double ek1, double ek2, double ek3, double delta){
//   return std::conj(G_AB_up_minus(ek1, ek2, ek3, delta));
// }
// cx_double G_BA_down_minus(double ek1, double ek2, double ek3, double delta){
//   return std::conj(G_AB_down_minus(ek1, ek2, ek3, delta));
// }
// cx_double G_BA_up_plus(double ek1, double ek2, double ek3, double delta){
//   return std::conj(G_AB_up_plus(ek1, ek2, ek3, delta));
// }
// cx_double G_BA_down_plus(double ek1, double ek2, double ek3, double delta){
//   return std::conj(G_AB_down_plus(ek1, ek2, ek3, delta));
// }

// double G_BB(double E, double delta, double spin){
//   return 0.5 * ( 1. + spin * delta / E );
// }
// double G_BB_minus(double ek1, double ek2, double ek3, double delta, double spin){
//   double E = ek3 - sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_BB( E, delta, spin );
// }
// double G_BB_plus(double ek1, double ek2, double ek3, double delta, double spin){
//   double E = ek3 + sqrt(delta*delta + ek1*ek1 + ek2*ek2);
//   return G_BB( E, delta, spin );
// }
// double G_BB_up_minus(double ek1, double ek2, double ek3, double delta){
//   return G_BB_minus( ek1, ek2, ek3, delta, 1.0 );
// }
// double G_BB_down_minus(double ek1, double ek2, double ek3, double delta){
//   return G_BB_minus( ek1, ek2, ek3, delta, - 1.0 );
// }
// double G_BB_up_plus(double ek1, double ek2, double ek3, double delta){
//   return G_BB_plus( ek1, ek2, ek3, delta, 1.0 );
// }
// double G_BB_down_plus(double ek1, double ek2, double ek3, double delta){
//   return G_BB_plus( ek1, ek2, ek3, delta, - 1.0 );
// }


double calc_bk_up_in(double e_free, double delta){
  double Ek = eigenenergy_HF_minus( e_free, delta );
  return  e_free / sqrt( ( Ek - delta ) * ( Ek - delta ) + e_free * e_free );
}

double calc_bk_up_out(double e_free, double delta){
  double Ek = eigenenergy_HF_plus( e_free, delta );
  
  /* Special treatment on the Fermi surface */
  double e_eps = 1e-12;
  if ( std::abs( e_free ) <= e_eps ) { return 1.; }
  else {
    return  e_free / sqrt( ( Ek - delta ) * ( Ek - delta ) + e_free * e_free );
  }
}

double calc_ak_up_in(double e_free, double delta){
  double bk = calc_bk_up_in( e_free, delta );
  return sqrt( 1. - bk * bk );
}

double calc_ak_up_out(double e_free, double delta){
  double bk = calc_bk_up_out( e_free, delta );
  return sqrt( 1. - bk * bk );
}

double calc_bk_down_in(double e_free, double delta){
  return calc_ak_up_in( e_free, delta );
}

double calc_bk_down_out(double e_free, double delta){
  return calc_ak_up_out( e_free, delta );
}

double calc_ak_down_in(double e_free, double delta){
  return calc_bk_up_in( e_free, delta );
}

double calc_ak_down_out(double e_free, double delta){
  return calc_bk_up_out( e_free, delta );
}

cx_double larger_eigenvalue(cx_double A, cx_double B, cx_double D){
  return 0.5 * ( A + D + sqrt( std::conj( A - D ) * ( A - D ) + 4. * std::conj(B) * B ) );
}

double wave_vector_in_BZ(double k){
  /* Returning - M_PI <= k < M_PI */
  double k_eps = 1e-12;
  while( k > M_PI - k_eps ) k -= 2. * M_PI;
  while( k < - M_PI - k_eps ) k += 2. * M_PI;
  if ( std::abs( k ) < k_eps ) k = 0;
  if ( std::abs( k + M_PI ) < k_eps ) k = - M_PI;
  return k;
}

void add_to_sus_mat(cx_double& A, cx_double& B, cx_double& D, double e_free, double e_free2, double delta, cx_double omega){
  double e_eps = 1e-12;
  double E1 = eigenenergy_HF_plus( e_free2, delta );
  double E2 = eigenenergy_HF_minus( e_free, delta );
  cx_double diff_E1 = E1 - E2 + omega;
  cx_double diff_E2 = E1 - E2 - omega;
  
  double ak_up_in = calc_ak_up_in( e_free, delta );
  double ak_q_up_out = calc_ak_up_out( e_free2, delta );
  double ak_down_in = calc_ak_down_in( e_free, delta );
  double ak_q_down_out = calc_ak_down_out( e_free2, delta );
  double bk_up_in = calc_bk_up_in( e_free, delta );
  double bk_q_up_out = calc_bk_up_out( e_free2, delta );
  double bk_down_in = calc_bk_down_in( e_free, delta );
  double bk_q_down_out = calc_bk_down_out( e_free2, delta );
  
  double factor = 1.;
  
  /* Taking into account a half of the contribution from the Fermi surface */
  if ( std::abs( e_free ) <= e_eps ) {
    factor = 0.5;
  }
   
  A += factor * ( ak_up_in * ak_up_in * ak_q_down_out * ak_q_down_out / diff_E1 + ak_down_in * ak_down_in * ak_q_up_out * ak_q_up_out / diff_E2 );
  B += factor * ( ak_up_in * bk_up_in * ak_q_down_out * bk_q_down_out / diff_E1 + ak_down_in * bk_down_in * ak_q_up_out * bk_q_up_out / diff_E2 );
  D += factor * ( bk_up_in * bk_up_in * bk_q_down_out * bk_q_down_out / diff_E1 + bk_down_in * bk_down_in * bk_q_up_out * bk_q_up_out / diff_E2 );
}

void add_to_sus_mat2(hoppings const& ts, double mu, cx_double& A, cx_double& B, cx_double& C, cx_double& D, double qx, double qy, double qz, double kx, double ky, double kz, double delta, cx_double omega, bool zz){
  double eps = 1e-12;
  
  double ek1 = ts.ek1(kx, ky, kz);
  double ek2 = ts.ek2(kx, ky, kz);
  double ek3 = ts.ek3(kx, ky, kz);
  double Ek = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  
  /* Checking if the eigenenergy is below the chemical potential. */
  double mu_free = 0;  /* Assume at half filling */
  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
  if ( e_free > mu_free + eps ) return;
    
  /* Prefactor */
  double factor = 1.;
  double factor2 = 1.;
  
  /* On the zone boundary */
  if ( std::abs(e_free - mu_free) < eps ) {
    factor = 0.5;
  }
  
  double diff_x = wave_vector_in_BZ( kx - qx );
  double diff_y = wave_vector_in_BZ( ky - qy );
  double diff_z = wave_vector_in_BZ( kz - qz );
  
  // double e_free2 = energy_free_electron( 1., mu_free, diff_x, diff_y );  /* ad-hoc */
  
  double ek_q1 = ts.ek1(diff_x, diff_y, diff_z);
  double ek_q2 = ts.ek2(diff_x, diff_y, diff_z);
  double ek_q3 = ts.ek3(diff_x, diff_y, diff_z);
  double Ek_q = eigenenergy_HF_plus(ek_q1, ek_q2, ek_q3, delta);
  
  cx_double diff_E1 = Ek_q - Ek + omega;
  cx_double diff_E2 = Ek_q - Ek - omega;
  
  cx_double ak_up = calc_ak_up_in_minus(ek1, ek2, ek3, delta);
  cx_double ak_down = calc_ak_down_in_minus(ek1, ek2, ek3, delta);
  cx_double bk_up = calc_bk_up_in_minus(ek1, ek2, ek3, delta);
  cx_double bk_down = calc_bk_down_in_minus(ek1, ek2, ek3, delta);

  cx_double ak_q_up, ak_q_down, bk_q_up, bk_q_down;
  ak_q_up = calc_ak_up_in_plus( ek_q1, ek_q2, ek_q3, delta);
  ak_q_down = calc_ak_down_in_plus( ek_q1, ek_q2, ek_q3, delta);
  bk_q_up = calc_bk_up_in_plus( ek_q1, ek_q2, ek_q3, delta);
  bk_q_down = calc_bk_down_in_plus( ek_q1, ek_q2, ek_q3, delta);

  // ak_q_up = calc_ak_up_in_plus( - ek_q1, - ek_q2, ek_q3, delta);
  // ak_q_down = calc_ak_down_in_plus( - ek_q1, - ek_q2, ek_q3, delta);
  // bk_q_up = calc_bk_up_in_plus( - ek_q1, - ek_q2, ek_q3, delta);
  // bk_q_down = calc_bk_down_in_plus( - ek_q1, - ek_q2, ek_q3, delta);
  
  if ( zz ) {
    /* Correct? */
    A += factor * factor2 * ( std::norm(ak_up) * std::norm(ak_q_up) / diff_E1 + std::norm(ak_up) * std::norm(ak_q_up) / diff_E2 );
    A += factor * factor2 * ( std::norm(ak_down) * std::norm(ak_q_down) / diff_E1 + std::norm(ak_down) * std::norm(ak_q_down) / diff_E2 );    
    
    B += factor * factor2 * ( std::conj(ak_up) * bk_up * ak_q_up * std::conj(bk_q_up) / diff_E1 + std::conj(ak_up) * bk_up * ak_q_up * std::conj(bk_q_up) / diff_E2 );
    B += factor * factor2 * ( std::conj(ak_down) * bk_down * ak_q_down * std::conj(bk_q_down) / diff_E1 + std::conj(ak_down) * bk_down * ak_q_down * std::conj(bk_q_down) / diff_E2 );    
    
    C += factor * factor2 * ( ak_up * std::conj(bk_up) * std::conj(ak_q_up) * bk_q_up / diff_E1 + ak_up * std::conj(bk_up) * std::conj(ak_q_up) * bk_q_up / diff_E2 );
    C += factor * factor2 * ( ak_down * std::conj(bk_down) * std::conj(ak_q_down) * bk_q_down / diff_E1 + ak_down * std::conj(bk_down) * std::conj(ak_q_down) * bk_q_down / diff_E2 );    

    D += factor * factor2 * ( std::norm(bk_up) * std::norm(bk_q_up) / diff_E1 + std::norm(bk_up) * std::norm(bk_q_up) / diff_E2 );
    D += factor * factor2 * ( std::norm(bk_down) * std::norm(bk_q_down) / diff_E1 + std::norm(bk_down) * std::norm(bk_q_down) / diff_E2 );    
  } else {
    A += factor * factor2 * ( std::norm(ak_up) * std::norm(ak_q_down) / diff_E1 + std::norm(ak_down) * std::norm(ak_q_up) / diff_E2 );
    
    B += factor * factor2 * ( std::conj(ak_up) * bk_up * ak_q_down * std::conj(bk_q_down) / diff_E1 + std::conj(ak_down) * bk_down * ak_q_up * std::conj(bk_q_up) / diff_E2 );
    
    C += factor * factor2 * ( ak_up * std::conj(bk_up) * std::conj(ak_q_down) * bk_q_down / diff_E1 + ak_down * std::conj(bk_down) * std::conj(ak_q_up) * bk_q_up / diff_E2 );

    D += factor * factor2 * ( std::norm(bk_up) * std::norm(bk_q_down) / diff_E1 + std::norm(bk_down) * std::norm(bk_q_up) / diff_E2 );
  }
}

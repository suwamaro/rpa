/*****************************************************************************
*
* Functions for mean field calculations.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "mean_field.h"
#include "rpa_util.h"

/* Member functions of MeanField */
MeanField::MeanField(std::size_t n):E_(n), U_(n,n){}


/* Member functions of MeanFieldBilayer */
MeanFieldBilayer::MeanFieldBilayer():MeanField(4){}

void MeanFieldBilayer::set_parameters(hoppings_bilayer2 const& h, double delta, double* kvec){
  hb_ = h;
  delta_ = delta;  
  for(int nu=0; nu < 3; nu++){ kvec_[nu] = kvec[nu]; }
}

void MeanFieldBilayer::set_eigen_energies(){
  double kx = kvec_[0];
  double ky = kvec_[1];
  double kz = kvec_[2];  
  cx_double ek1 = hb_.ek1(kx, ky, kz);
  cx_double ek23 = hb_.ek23(kx, ky, kz);
  cx_double ekz = hb_.ekz(kx, ky, kz);
  double ek_up = eigenenergy_HF(up_spin, ek1, ek23, ekz, hb_.tz, kz, delta_);
  double ek_down = eigenenergy_HF(down_spin, ek1, ek23, ekz, hb_.tz, kz, delta_);
  E_[0] = ek_up;
  E_[1] = ek_up;
  E_[2] = ek_down;
  E_[3] = ek_down;
}

void MeanFieldBilayer::set_eigen_vectors(){
  double kx = kvec_[0];
  double ky = kvec_[1];
  double kz = kvec_[2];    
  cx_double ek1 = hb_.ek1(kx, ky, kz);
  cx_double ek23 = hb_.ek23(kx, ky, kz);
  cx_double ekz = hb_.ekz(kx, ky, kz);    
  cx_double xki_up = xk(up_spin, ek1, hb_.tz, kz, delta_);
  cx_double xki_down = xk(down_spin, ek1, hb_.tz, kz, delta_);  
  double zki = zk(ek1, hb_.tz, kz, delta_);

  U_.zeros();
  U_.col(0) = calc_eigen_vector(1, up_spin, xki_up, zki);
  U_.col(1) = calc_eigen_vector(1, down_spin, xki_down, zki);
  U_.col(2) = calc_eigen_vector(-1, up_spin, xki_up, zki);
  U_.col(3) = calc_eigen_vector(-1, down_spin, xki_down, zki);
  // U_dagger_ = arma::conj(U_);
}

cx_vec MeanFieldBilayer::calc_eigen_vector(int sign, int spin, cx_double xk, double zk) const {
  cx_vec X(4, arma::fill::zeros);
  if ( spin == up_spin ) {
    X[0] = xk * sqrt(0.5 * (1. + (double)sign * zk));
    X[2] = (double)sign * sqrt(0.5 * (1. - (double)sign * zk));
  } else /* spin == down_spin */ {
    X[1] = xk * sqrt(0.5 * (1. - (double)sign * zk));
    X[3] = (double)sign * sqrt(0.5 * (1. + (double)sign * zk));
  }
  return X;
}

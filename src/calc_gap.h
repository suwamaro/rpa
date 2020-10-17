/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __CALC_GAP__
#define __CALC_GAP__

#include <armadillo>
#include "rpa.h"
#include "parameters.h"
#include "cuba_helper.h"

/* For the bare response function */
cx_double calc_prefactor_bare_res_func_bilayer(int sg1, int sg2, hoppings2 const& ts, double kx, double ky, double kz, double qx, double qy, double qz, cx_double omega, double delta);

/* Coefficients */
cx_double calc_ak_up_in_minus(double ek1, double ek2, double ek3, double delta);
cx_double calc_ak_down_in_minus(double ek1, double ek2, double ek3, double delta);

/* Common to the lattices */
cx_double larger_eigenvalue(cx_double A, cx_double B, cx_double D);
double wave_vector_in_BZ(double k);
void add_to_sus_mat(cx_double& A, cx_double& B, cx_double& D, double e_free, double e_free2, double delta, cx_double omega);
void add_to_sus_mat2(hoppings const& ts, double mu, cx_double& A, cx_double& B, cx_double& C, cx_double& D, double qx, double qy, double qz, double kx, double ky, double kz, double delta, cx_double omega, bool zz);
void add_to_sus_mat3(hoppings2 const& ts, double mu, arma::cx_mat& chi, double qx, double qy, double qz, double kx, double ky, double kz, double delta, cx_double omega, bool zz);

struct Polarization {
  Polarization(){};
  explicit Polarization(int Lx, int Ly, int Lz, int nsub);
  void set_q(double qx, double qy, double qz);
  void set_table(hoppings2 const& ts, double delta);
  int nbands() const { return nbands_; }
  bool is_table_set() const { return is_table_set_; }
  int Lx() const { return Lx_; }
  int Ly() const { return Ly_; }
  int Lz() const { return Lz_; }  
  int nsub() const { return nsub_; }
  double qx() const { return qx_; }
  double qy() const { return qy_; }
  double qz() const { return qz_; }
  long unsigned int table_size() const;
  
  /* Assume that k = 2pi / L * m, where m is an integer. */
  int pullback(double k, int L) const;
  int xyz_to_index(int x, int y, int z) const;
  int k_to_index(double kx, double ky, double kz) const;
  void calc_polarization(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* ppm, cx_double* pzz) const;  
  void build_table(hoppings2 const& ts, double delta);
  void get_Ppm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Ppmk) const;
  void get_Pzz(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Pzzk) const;  
  ~Polarization();

  int nbands_ = NSUBL;
  bool is_table_set_;
  cx_double *Ppm_ = nullptr;
  cx_double *Pzz_ = nullptr;
  int Lx_;
  int Ly_;
  int Lz_;
  int nsub_;
  double qx_;
  double qy_;
  double qz_;
};

void add_to_sus_mat4(hoppings2 const& ts, double mu, arma::cx_mat& chi_pm, arma::cx_mat& chi_zz, double kx, double ky, double kz, Polarization const& pz, double delta, cx_double omega);

/* For a square lattice */
double calc_eigval_square(int L, double t, double mu, double U, double delta, double qx, double qy, double omega);
double calc_gap_square(int L, double t, double mu, double U, double delta, double qx, double qy);
cx_double calc_intensity_square(int L, hoppings const& ts, double mu, double U, double delta, double qx, double qy, cx_double omega, bool zz);

/* For a bilayer lattice */
cx_double calc_intensity_bilayer(int L, hoppings const& ts, double mu, double U, double delta, double qx, double qy, double qz, cx_double omega, bool zz);
std::tuple<cx_double, cx_double> calc_intensity_bilayer2(int L, hoppings_bilayer2& ts, double mu, double U, double delta, CubaParam const& cbp, Polarization const& pz, cx_double omega, bool continuous_k);

/* For a simple cubic lattice */
cx_double calc_intensity_cubic(int L, hoppings const& ts, double mu, double U, double delta, double qx, double qy, double qz, cx_double omega, bool zz);
double calc_eigval_cubic(int L, double t, double mu, double U, double delta, double qx, double qy, double qz, double omega);
double calc_gap_cubic(int L, double t, double mu, double U, double delta, double qx, double qy, double qz);

# endif // __CALC_GAP__

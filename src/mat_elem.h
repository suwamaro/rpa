/*****************************************************************************
*
* Functions for matrix elements
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _MAT_ELEM_H
#define _MAT_ELEM_H

#include "rpa.h"
#include "array3.h"

typedef rpa::array3 vec3;

class MatElem {
public:
  MatElem(){};
  explicit MatElem(int Lx, int Ly, int Lz);
  void set_q(double qx, double qy, double qz);
  void assign_is_table_set(bool b) { is_table_set_ = b; }
  int nbands() const { return nbands_; }
  bool is_table_set() const { return is_table_set_; }
  int Lx() const { return Lx_; }
  int Ly() const { return Ly_; }
  int Lz() const { return Lz_; }
  double qx() const { return qx_; }
  double qy() const { return qy_; }
  double qz() const { return qz_; }
  
  /* Assume that k = 2pi / L * m, where m is an integer. */
  int pullback(double k, int L) const;
  int xyz_to_index(int x, int y, int z) const;
  int k_to_index(double kx, double ky, double kz) const;
  void build_table(hoppings2 const& ts, double delta);
  ~MatElem(){};

private:
  // OperatorK opek;
  int nbands_ = NSUBL;
  bool is_table_set_;
  int Lx_;
  int Ly_;
  int Lz_;
  double qx_;
  double qy_;
  double qz_;
};

class BasicMatrix {
public:
  BasicMatrix();
  mat tau_A, tau_B;
  cx_mat Pauli_0, Pauli_p, Pauli_m, Pauli_z;  
  cx_mat Sigma_A0, Sigma_Ap, Sigma_Am, Sigma_Az;
  cx_mat Sigma_B0, Sigma_Bp, Sigma_Bm, Sigma_Bz;  
};

/* Matrix elements for the charge and magnetic susceptibility */
class MatElemF : public MatElem, public BasicMatrix {
public:
  MatElemF(){};
  explicit MatElemF(int Lx, int Ly, int Lz, int nsub);
  int nsub() const { return nsub_; }  
  void set_table(hoppings2 const& ts, double delta);
  std::size_t table_size() const;  
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* F00_up, cx_double* F00_down, cx_double* Fpm, cx_double* Fzz_up, cx_double* Fzz_down) const;  
  void build_table(hoppings2 const& ts, double delta);
  void get_00_up(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const;
  void get_00_down(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const;  
  void get_pm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fpmk) const;
  void get_zz_up(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const;
  void get_zz_down(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const;    
  ~MatElemF();

private:
  int nsub_;  
  // OperatorK opek;
  cx_double *F00_up_ = nullptr;
  cx_double *F00_down_ = nullptr;  
  cx_double *Fpm_ = nullptr;
  cx_double *Fzz_up_ = nullptr;
  cx_double *Fzz_down_ = nullptr;  
};

/* Matrix elements for the effective Raman operator */
class MatElemK : public MatElem {
public:
  MatElemK(){};
  explicit MatElemK(int Lx, int Ly, int Lz, int mu, int nu, std::vector<BondDelta> const& bonds);
  void set_table(hoppings2 const& ts, double delta);
  std::size_t table_size() const;
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double* R) const;  
  void build_table(hoppings2 const& ts, double delta);
  void get_elem(double kx, double ky, double kz, cx_double* R) const;
  ~MatElemK();

private:
  BondDelta mu_;
  BondDelta nu_;  
  std::vector<BondDelta> bonds_;
  cx_double *R_up_ = nullptr;
};

#endif // _MAT_ELEM_H

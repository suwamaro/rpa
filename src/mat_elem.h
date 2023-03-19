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
  explicit MatElem(int Lx, int Ly, int Lz, int u_cell);
  explicit MatElem(int Lx, int Ly, int Lz, int u_cell, int n_coefs);  
  void set_q(double qx, double qy, double qz);
  void assign_is_table_set(bool b) { is_table_set_ = b; }
  int n_bands() const { return n_bands_; }
  int n_spins() const { return n_spins_; }  
  int n_coefs() const { return n_coefs_; }
  bool is_table_set() const { return is_table_set_; }
  std::size_t system_size() const { return (std::size_t)Lx() * Ly() * Lz(); }
  int Lx() const { return Lx_; }
  int Ly() const { return Ly_; }
  int Lz() const { return Lz_; }
  int u_cell_size() const { return u_cell_size_; }  
  double qx() const { return qx_; }
  double qy() const { return qy_; }
  double qz() const { return qz_; }
  
  /* Assume that k = 2pi / L * m, where m is an integer. */
  int pullback(double k, int L) const;
  int xyz_to_index(int x, int y, int z) const;
  std::size_t k_to_index(double kx, double ky, double kz) const;

  /* Virtual functions */
  virtual void set_table(hoppings2 const& ts, double delta) = 0;
  virtual std::size_t table_size() const = 0;
  virtual void build_table(hoppings2 const& ts, double delta) = 0;
  virtual void get_elem(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *M) const {};
  virtual ~MatElem(){};

private:
  // OperatorK opek;
  int n_bands_ = NSUBL;
  int n_spins_ = 2;  
  int n_coefs_;  
  bool is_table_set_;
  int Lx_;
  int Ly_;
  int Lz_;
  int u_cell_size_;
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
  virtual ~BasicMatrix(){};
};

/* Matrix elements for the charge and magnetic susceptibility */
class MatElemF : public MatElem, public BasicMatrix {
public:
  MatElemF(){};
  explicit MatElemF(int Lx, int Ly, int Lz, int u_cell, int n_sub);
  int n_sub() const { return n_sub_; }  
  void set_table(hoppings2 const& ts, double delta) override;
  std::size_t table_size() const override;  
  void build_table(hoppings2 const& ts, double delta) override;
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* F00_up, cx_double* F00_down, cx_double* Fpm, cx_double* Fzz_up, cx_double* Fzz_down) const;  
  void get_00_up(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const;
  void get_00_down(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const;  
  void get_pm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fpmk) const;
  void get_zz_up(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const;
  void get_zz_down(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const;    
  ~MatElemF();

private:
  int n_sub_;  
  // OperatorK opek;
  cx_double *F00_up_ = nullptr;
  cx_double *F00_down_ = nullptr;  
  cx_double *Fpm_ = nullptr;
  cx_double *Fzz_up_ = nullptr;
  cx_double *Fzz_down_ = nullptr;  
};

/* Matrix elements for the current operator */
class MatElemVelocity : public MatElem {
public:
  MatElemVelocity(){};
  explicit MatElemVelocity(int Lx, int Ly, int Lz, int u_cell, std::vector<BondDelta> const& bonds);
  int n_directions() const { return n_directions_; }
  void set_table(hoppings2 const& ts, double delta) override;
  std::size_t table_size() const override;
  void build_table(hoppings2 const& ts, double delta) override;
  int direction_index(BondDelta const& mu) const;
  void calc_mat_elems(hoppings2 const& ts, cx_mat const& Udg, cx_mat const& U, cx_mat const& U_bar, double kx, double ky, double kz, BondDelta const& e_mu, int sign1, int sign2, int sigma, cx_double *M) const;
  std::size_t get_index(std::size_t xyz, int bands, int dir, int spin) const;
  int bands_index(int sign1, int sign2) const;
  int spin_index(int sigma) const;
  void get_elem(double kx, double ky, double kz, int sign1, int sign2, int dir, int spin, cx_double *M) const;  
  ~MatElemVelocity();

private:
  int n_directions_ = 3;
  vec3 q0_;  
  std::vector<BondDelta> bonds_;
  cx_double *velocity_ = nullptr;
};

/* Matrix elements for the resonant contributions of the effective Raman operator */
class MatElemK : public MatElem {
public:
  MatElemK(){};
  explicit MatElemK(int Lx, int Ly, int Lz, int u_cell, int mu, int nu, std::vector<BondDelta> const& bonds);
  void set_table(hoppings2 const& ts, double delta) override;
  void set_table(hoppings2 const& ts, double delta, MatElemVelocity const& mev);  
  std::size_t table_size() const override;
  void build_table(hoppings2 const& ts, double delta) override;
  void build_table(hoppings2 const& ts, double delta, MatElemVelocity const& mev);  
  void get_elem(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *R) const override;
  void set_occupied_and_empty_vectors(hoppings2 const& ts, double delta, double ch_pot);  
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *R) const;
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, MatElemVelocity const& mev, cx_double *R) const;  
  ~MatElemK();

private:
  BondDelta mu_;
  BondDelta nu_;
  std::vector<BondDelta> bonds_;
  cx_double *R_up_ = nullptr;
  std::vector<vec3> occupied_, empty_;   // Occupied and empty wave vectors
  std::vector<std::size_t> not_half_occupied_;   // The indices of the occupied and empty vectors.
};

/* Matrix elements for the nonresonant contributions of the effective Raman operator */
class MatElemN : public MatElem {
public:
  MatElemN(){};
  explicit MatElemN(int Lx, int Ly, int Lz, int u_cell, int mu, int nu, std::vector<BondDelta> const& bonds);  
  void set_table(hoppings2 const& ts, double delta) override;
  std::size_t table_size() const override;
  void build_table(hoppings2 const& ts, double delta) override;
  void get_elem(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *N) const override;
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *N) const;  
  ~MatElemN();

private:
  BondDelta mu_;
  BondDelta nu_;  
  std::vector<BondDelta> bonds_;
  cx_double *N_ = nullptr;
};

#endif // _MAT_ELEM_H

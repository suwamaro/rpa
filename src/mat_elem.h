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

struct MatElem {
  MatElem(){};
  explicit MatElem(int Lx, int Ly, int Lz, int nsub);
  void set_q(double qx, double qy, double qz);
  int nbands() const { return nbands_; }
  bool is_table_set() const { return is_table_set_; }
  int Lx() const { return Lx_; }
  int Ly() const { return Ly_; }
  int Lz() const { return Lz_; }  
  int nsub() const { return nsub_; }
  double qx() const { return qx_; }
  double qy() const { return qy_; }
  double qz() const { return qz_; }
  
  /* Assume that k = 2pi / L * m, where m is an integer. */
  int pullback(double k, int L) const;
  int xyz_to_index(int x, int y, int z) const;
  int k_to_index(double kx, double ky, double kz) const;
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* p00, cx_double* ppm, cx_double* pzz) const;  
  void build_table(hoppings2 const& ts, double delta);
  ~MatElem(){};

  // OperatorK opek;
  int nbands_ = NSUBL;
  bool is_table_set_;
  int Lx_;
  int Ly_;
  int Lz_;
  int nsub_;
  double qx_;
  double qy_;
  double qz_;
};

struct MatElemF : public MatElem {
  MatElemF(){};
  explicit MatElemF(int Lx, int Ly, int Lz, int nsub);
  void set_table(hoppings2 const& ts, double delta);
  std::size_t table_size() const;
  
  /* Assume that k = 2pi / L * m, where m is an integer. */
  void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* p00, cx_double* ppm, cx_double* pzz) const;  
  void build_table(hoppings2 const& ts, double delta);
  void get_00(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const;  
  void get_pm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fpmk) const;
  void get_zz(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const;  
  ~MatElemF();

  // OperatorK opek;
  cx_double *F00_ = nullptr;  
  cx_double *Fpm_ = nullptr;
  cx_double *Fzz_ = nullptr;
};

// struct MatElemK {
//   MatElemK(){};
//   explicit MatElemK(int Lx, int Ly, int Lz, int nsub);
//   void set_q(double qx, double qy, double qz);
//   void set_table(hoppings2 const& ts, double delta);
//   int nbands() const { return nbands_; }
//   bool is_table_set() const { return is_table_set_; }
//   int Lx() const { return Lx_; }
//   int Ly() const { return Ly_; }
//   int Lz() const { return Lz_; }  
//   int nsub() const { return nsub_; }
//   double qx() const { return qx_; }
//   double qy() const { return qy_; }
//   double qz() const { return qz_; }
//   std::size_t table_size() const;
  
//   /* Assume that k = 2pi / L * m, where m is an integer. */
//   int pullback(double k, int L) const;
//   int xyz_to_index(int x, int y, int z) const;
//   int k_to_index(double kx, double ky, double kz) const;
//   void calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* p00, cx_double* ppm, cx_double* pzz) const;  
//   void build_table(hoppings2 const& ts, double delta);
//   void get_00(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* K00k) const;  
//   void get_pm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Kpmk) const;
//   void get_zz(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Kzzk) const;  
//   ~MatElemK();

//   // OperatorK opek;
//   int nbands_ = NSUBL;
//   bool is_table_set_;
//   cx_double *K00_ = nullptr;  
//   cx_double *Kpm_ = nullptr;
//   cx_double *Kzz_ = nullptr;
//   int Lx_;
//   int Ly_;
//   int Lz_;
//   int nsub_;
//   double qx_;
//   double qy_;
//   double qz_;
// };

#endif // _MAT_ELEM_H

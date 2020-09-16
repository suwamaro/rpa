/*****************************************************************************
*
* RPA calculation for the fermionic Hubbard model
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __HOPPINGS__
#define __HOPPINGS__

#include "rpa.h"

class hoppings {
public:
  hoppings(){};  
  virtual double ek1(double kx, double ky, double kz) const = 0;
  virtual double ek2(double kx, double ky, double kz) const = 0;
  virtual double ek3(double kx, double ky, double kz) const = 0;
  virtual ~hoppings(){};
  double t_max() const;
  
  /* Parameters of an effective one-band model */
  double t = 0;
  double t_bar = 0;
  double tp = 0;
  double tpp = 0;
  double tz = 0;
  double tz_bar = 0;
  double tzp = 0;
};

class hoppings2 {
public:
  hoppings2(){};
  virtual cx_double ek1(double kx, double ky, double kz) const = 0;
  virtual cx_double ek23(double kx, double ky, double kz) const = 0;
  virtual cx_double ekz(double kx, double ky, double kz) const = 0;
  virtual ~hoppings2(){};
  double t_max() const;
  
  /* Parameters of an effective one-band model */
  cx_double t = 0;
  cx_double tp = 0;
  cx_double tpp = 0;
  cx_double tz = 0;
  cx_double tzp = 0;
};

class hoppings_square : public hoppings {
public:
  double ek1(double kx, double ky, double kz = 0) const;
  double ek2(double kx, double ky, double kz = 0) const;
  double ek3(double kx, double ky, double kz = 0) const;
  explicit hoppings_square(double v, double v_bar);

  // instantiations
  static std::unique_ptr<hoppings_square> mk_square(double v, double v_bar);
};

class hoppings_cubic : public hoppings {
public:
  double ek1(double kx, double ky, double kz) const;
  double ek2(double kx, double ky, double kz) const;
  double ek3(double kx, double ky, double kz) const;
  explicit hoppings_cubic(double t);

  // instantiations
  static std::unique_ptr<hoppings_cubic> mk_cubic(double v);
};

class hoppings_bilayer : public hoppings {
public:
  hoppings_bilayer(){};
  explicit hoppings_bilayer(double v, double v_bar, double vp, double vpp, double vz, double vz_bar, double vzp);  
  double ek1(double kx, double ky, double kz) const;
  double ek2(double kx, double ky, double kz) const;
  double ek3(double kx, double ky, double kz) const;

  // instantiations
  static std::unique_ptr<hoppings_bilayer> mk_bilayer(double v, double v_bar, double vp, double vpp, double vz, double vz_bar, double vzp);
};

class hoppings_bilayer2 : public hoppings2 {
public:
  hoppings_bilayer2(){};
  virtual ~hoppings_bilayer2(){}
  explicit hoppings_bilayer2(cx_double v, cx_double vp, cx_double vpp, cx_double vz, cx_double vzp);
  cx_double ek1(double kx, double ky, double) const;
  cx_double ek23(double kx, double ky, double) const;
  cx_double ekz(double kx, double ky, double kz) const;

  // instantiations
  static std::unique_ptr<hoppings_bilayer2> mk_bilayer2(cx_double v, cx_double vp, cx_double vpp, cx_double vz, cx_double vzp);
};

class hoppings_Sr3Ir2O7 : public hoppings_bilayer {
public:  
  explicit hoppings_Sr3Ir2O7(double theta, double phi, double t3);
  
  // instantiations
  static std::unique_ptr<hoppings_Sr3Ir2O7> mk_Sr3Ir2O7(double theta, double phi, double t3);
  
private:
  /* Parameters of the t2g three-band model */
  /* Reference: J.-M. Carter, et al., PRB 87 014433 (2013) 
                S. Mohapatra, et al., PRB 95, 094435 (2017). */  
  double t1 = 0.2;
  double t2 = -0.032;
  double t3 = -0.03;     // Modified
  // double t3 = -0.02  // in Reference
  double t4 = 0.188;
  double t5 = -0.054;
  double t6 = 0;
  double tm = -0.03;
  double tmp = 0.022;
  double t1z = 0.03;
  double t4z = -0.16;
  double tmz = 0.072;
  double t6z = -0.04;
  double t2z = 0;
};

#endif // __HOPPINGS__

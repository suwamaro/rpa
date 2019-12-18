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

class hoppings {
public:
  double ek1(double kx, double ky, double kz) const;
  double ek2(double kx, double ky, double kz) const;
  double ek3(double kx, double ky, double kz) const;
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

class hoppings_simple : public hoppings {
public:
  hoppings_simple(double t);
};

class hoppings_Sr3Ir2O7 : public hoppings {
public:  
  hoppings_Sr3Ir2O7(double theta, double phi, double t3);

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

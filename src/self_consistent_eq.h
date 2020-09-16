/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"

/* Integrand */
class SelfConsistentIntegrand {
public:
  SelfConsistentIntegrand(){}
  void set_parameters(hoppings2 *ts, double delta);
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  hoppings2 *ts() const { return ts_; }
  double delta() const { return delta_; }
  
private:  
  hoppings2 *ts_;
  double delta_;
};

class SelfConsistentIntegrandBilayer : public SelfConsistentIntegrand {
public:
  SelfConsistentIntegrandBilayer():SelfConsistentIntegrand(){}  
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
};

/* For a square lattice */
double self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double delta);
double solve_self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double U);

/* For a bilayer lattice */
double self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double delta);
double solve_self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double U);
double solve_self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double U, CubaParam const& cbp, bool continuous_k);

/* For a simple cubic lattice */
double self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double delta);
double solve_self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double U);

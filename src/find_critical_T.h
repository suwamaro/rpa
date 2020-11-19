/*****************************************************************************
*
* Functions for calculating the critical T.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"

/* Integrand */
class FindTcIntegrand {
public:
  explicit FindTcIntegrand(hoppings2 *ts):ts_(ts){}
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  void set_T(double _T) { T_ = _T; }
  hoppings2 *ts() const { return ts_; }
  double T() const { return T_; }
  
private:  
  hoppings2 *ts_;
  double T_;
};

class FindTcIntegrandBilayer : public FindTcIntegrand {
public:
  FindTcIntegrandBilayer():FindTcIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& h, double T);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;

private:
  hoppings_bilayer2 hb_;
};

/* For bilayer lattices */
void find_critical_T_bilayer(path& base_dir, rpa::parameters const& pr);

/*****************************************************************************
*
* Functions for calculating the dispersion
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"

/* Integrand */
class FindUcIntegrand {
public:
  explicit FindUcIntegrand(hoppings2 *ts):ts_(ts){}
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  hoppings2 *ts() const { return ts_; }
  
private:  
  hoppings2 *ts_;  
};

class FindUcIntegrandBilayer : public FindUcIntegrand {
public:
  FindUcIntegrandBilayer():FindUcIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& h);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;

private:
  hoppings_bilayer2 hb_;
};

/* For bilayer lattices */
void find_critical_U_bilayer(rpa::parameters const& pr);
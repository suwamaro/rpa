/*****************************************************************************
*
* Functions for finding the critical point.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "hoppings.h"
#include "cuba_helper.h"

/* Integrand */
class FindQCPIntegrand {
public:
  explicit FindQCPIntegrand(hoppings2 *ts):ts_(ts){}
  virtual double integrand(const double *qvec) const = 0;
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  void set_mu(double _mu);
  hoppings2 *ts() const;
  double mu() const;
  
private:  
  hoppings2 *ts_;
  double mu_;
};

class FindQCPIntegrandBilayer : public FindQCPIntegrand {
public:
  FindQCPIntegrandBilayer():FindQCPIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& h, double _mu);
  double integrand(const double *qvec) const;  
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;

private:
  hoppings_bilayer2 hb_;
};

/* For bilayer lattices */
double find_critical_point_bilayer(rpa::parameters& pr);
void find_critical_point_bilayer_output(path& base_dir, rpa::parameters& pr);

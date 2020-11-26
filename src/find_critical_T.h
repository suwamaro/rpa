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
  virtual double integrand(const double *qvec) const = 0;
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  void set_T(double _T);
  void set_mu(double _mu);
  hoppings2 *ts() const;
  double T() const;
  double mu() const;
  
private:  
  hoppings2 *ts_;
  double T_;
  double mu_;
};

class FindTcIntegrandBilayer : public FindTcIntegrand {
public:
  FindTcIntegrandBilayer():FindTcIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& h, double _T, double _mu);
  double integrand(const double *qvec) const;  
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;

private:
  hoppings_bilayer2 hb_;
};

/* For bilayer lattices */
double find_critical_T_bilayer(rpa::parameters const& pr);
void find_critical_T_bilayer_output(path& base_dir, rpa::parameters const& pr);

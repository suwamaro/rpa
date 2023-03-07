/*****************************************************************************
*
* Functions for calculating the critical U.
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
class FindUcIntegrand {
public:
  explicit FindUcIntegrand(hoppings2 *ts):ts_(ts){}
  virtual double integrand(const double *qvec) const = 0;
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  hoppings2 *ts() const { return ts_; }
  
private:  
  hoppings2 *ts_;
};

class FindUcIntegrandBilayer : public FindUcIntegrand {
public:
  FindUcIntegrandBilayer():FindUcIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& h, double delta, double mu);
  double integrand(const double *qvec) const;  
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
  double delta() const { return delta_; }  
  double mu() const { return mu_; }
  
private:
  hoppings_bilayer2 hb_;
  double delta_;
  double mu_;  
};

/* For bilayer lattices */
double find_critical_U_bilayer(rpa::parameters const& pr);
void find_critical_U_bilayer_output(path& base_dir, rpa::parameters const& pr);
void check_mean_field_eq_bilayer(std::string const& ofilen, rpa::parameters const& pr);
std::tuple<double, double, double, double, double, double> calc_total_energies(rpa::parameters const& pr, double U, double delta_i, double mu_i, bool set_initial_values);

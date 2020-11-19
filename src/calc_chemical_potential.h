/*****************************************************************************
*
* Chemical potential calculation for the fermionic Hubbard model.
*
* Copyright (C) 2019 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __CALC_CHEMICAL_POTENTIAL__
#define __CALC_CHEMICAL_POTENTIAL__

#include <cuba.h>
#include "cuba_helper.h"
#include "rpa_util.h"

/* Integrand */
class ElecFillingIntegrand {
public:
  explicit ElecFillingIntegrand(hoppings2 *ts):ts_(ts){}
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  virtual void set_parameters(double _mu, double _T, double _delta);
  hoppings2 *ts() const { return ts_; }
  double mu() const { return mu_; }
  double T() const { return T_; }
  double delta() const { return delta_; }  
  
private:  
  hoppings2 *ts_;
  double mu_;
  double T_;
  double delta_;  
};

class ElecFillingIntegrandBilayer : public ElecFillingIntegrand {
public:
  ElecFillingIntegrandBilayer():ElecFillingIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& ts, double mu, double T, double delta);  
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
  
private:
  hoppings_bilayer2 hb_;
};

double calc_chemical_potential_bilayer(int L, hoppings const& ts, double delta);
double calc_chemical_potential_bilayer2(path& base_dir, int L, hoppings2 const& ts, double delta);
double calc_chemical_potential_bilayer3(path& base_dir, int L, hoppings_bilayer2 const& ts, double filling, double T, double delta, CubaParam const& cbp, bool continuous_k);
std::tuple<double, double> calc_charge_gap_bilayer(int L, hoppings2 const& ts, double delta);

#endif // __CALC_CHEMICAL_POTENTIAL__

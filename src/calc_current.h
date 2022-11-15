/*****************************************************************************
*
* Functions for calculating electronic current.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"
#include "mean_field.h"

/* Integrand */
class CurrentIntegrand {
public:
  explicit CurrentIntegrand(hoppings2 *ts):ts_(ts){}
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) = 0;
  hoppings2 *ts() const { return ts_; }
  constexpr int x_dir() const { return 0; }
  constexpr int y_dir() const { return 1; }
  constexpr int z_dir() const { return 2; }  
  virtual ~CurrentIntegrand(){}
  
private:  
  hoppings2 *ts_;
};

class CurrentIntegrandBilayer : public CurrentIntegrand {
public:
  CurrentIntegrandBilayer():CurrentIntegrand(&hb_),mfb_(){}
  void set_parameters(hoppings_bilayer2 const& h, double kT, double mu_, double delta, BondDelta bd, int g1, int g2, double *Qvec);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);
  void set_variables(double *ks);
  cx_double integrand(double *ks);
  double kT() const { return kT_; }
  double mu() const { return mu_; }
  double delta() const { return delta_; }
  int gamma1() const { return gamma1_; }
  int gamma2() const { return gamma2_; }
  cx_double t() const { return t_; }
  cx_double phase() const { return phase_; }
  cx_double epsilon(int m) const { return epsilon_[m]; }
  bool no_contribution() const { return std::abs(t()) < 1e-12; };
  
private:
  hoppings_bilayer2 hb_;
  MeanFieldBilayer mfb_;
  double kT_;
  double mu_;
  double delta_;
  BondDelta bond_;
  int gamma1_, gamma2_;    
  double Qvec_[3];
  cx_double t_;
  cx_double phase_;
  cx_double epsilon_[2];  // Sublattice factor for the two terms.
  cx_mat UfUd_;
};

/* For bilayer lattices */
void calc_current_bilayer(path& base_dir, rpa::parameters const& pr);

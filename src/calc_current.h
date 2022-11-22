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
  cx_double hopping_U1(int g1, int g2, int sigma1, int sigma2);
  cx_double integrand(double *ks);
  double kT() const { return kT_; }
  double mu() const { return mu_; }
  double delta() const { return delta_; }
  int gamma1() const { return gamma1_; }
  int gamma2() const { return gamma2_; }
  cx_double epsilon() const { return epsilon_; }
  bool no_contribution() const;
  
private:
  hoppings_bilayer2 hb_;
  MeanFieldBilayer mfb_;
  double kT_;
  double mu_;
  double delta_;
  BondDelta bond_;
  int gamma1_, gamma2_;    
  double Qvec_[3];
  cx_double epsilon_;  // Sublattice factor
  cx_mat UfUd_;
};

/* For bilayer lattices */
void calc_current_bilayer(path& base_dir, rpa::parameters const& pr);

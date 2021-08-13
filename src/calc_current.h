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
  CurrentIntegrandBilayer():CurrentIntegrand(&hb_),mfb_(), J_(4,4){}
  void set_parameters(hoppings_bilayer2 const& h, double kT, double mu_, double delta, int dir_, double *Qvec);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);
  void set_current_matrix(double *kvec, int dir);
  cx_double integrand(double kx, double ky, double kz);
  double kT() const { return kT_; }
  double mu() const { return mu_; }
  double delta() const { return delta_; }
  int dir() const { return dir_; }

private:
  hoppings_bilayer2 hb_;
  MeanFieldBilayer mfb_;
  double kT_;
  double mu_;
  double delta_;
  int dir_;
  double Qvec_[3];
  cx_mat J_;
};

/* For bilayer lattices */
void calc_current_bilayer(path& base_dir, rpa::parameters const& pr);

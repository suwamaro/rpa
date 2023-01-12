/*****************************************************************************
*
* Functions for calculating the response function
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"
#include "calc_gap.h"
#include "mat_elem.h"

class ResponseFuncIntegrand {
public:
  explicit ResponseFuncIntegrand(hoppings2 *ts):ts_(ts){}
  virtual ~ResponseFuncIntegrand(){}
  void update_parameters(double T, double delta, double mu, cx_double omega, MatElemF const& me_F);
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;  
  hoppings2 *ts() const { return ts_; }
  double T() const { return T_; }
  double delta() const { return delta_; }
  double mu() const { return mu_; }
  cx_double omega() const { return omega_; }  
  const MatElemF *me_F() const { return &me_F_; }
  
private:  
  hoppings2 *ts_;
  double T_;
  double delta_;
  double mu_;
  cx_double omega_;
  MatElemF me_F_;  
};

class ResponseFuncIntegrandBilayer : public ResponseFuncIntegrand {
public:
  ResponseFuncIntegrandBilayer():ResponseFuncIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& h, double T, double delta, double mu, cx_double omega, MatElemF const& me_F);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;

private:
  hoppings_bilayer2 hb_;
};

std::tuple<arma::cx_mat, arma::cx_mat> calc_bare_response_bilayer(int L, hoppings_bilayer2 const& ts, double mu, double U, double T, double delta, CubaParam const& cbp, MatElemF const& me_F, cx_double omega, bool continuous_k);

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

class ResponseFuncIntegrand {
public:
  explicit ResponseFuncIntegrand(hoppings2 *ts):ts_(ts){}
  virtual ~ResponseFuncIntegrand(){}
  void update_parameters(double delta, cx_double omega, Polarization const& Pz);
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;  
  hoppings2 *ts() const { return ts_; }
  double delta() const { return delta_; }
  cx_double omega() const { return omega_; }  
  const Polarization *Pz() const { return &Pz_; }
  
private:  
  hoppings2 *ts_;
  double delta_;
  cx_double omega_;
  Polarization Pz_;  
};

class ResponseFuncIntegrandBilayer : public ResponseFuncIntegrand {
public:
  ResponseFuncIntegrandBilayer():ResponseFuncIntegrand(&hb_){}
  void set_parameters(hoppings_bilayer2 const& h, double delta, cx_double omega, Polarization const& Pz);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;

private:
  hoppings_bilayer2 hb_;
};

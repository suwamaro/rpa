/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "self_consistent_eq.h"
#include "rpa_util.h"

/* Member functions of SelfConsistentIntegrand */
void SelfConsistentIntegrand::set_parameters(hoppings2 *ts, double delta){
  ts_ = ts;
  delta_ = delta;
}

/* Member functions of SelfConsistentIntegrandBilayer */
int SelfConsistentIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
  /* Reset */
  ff[0] = 0;
  
  /* Wavenumbers */
  double k1 = xx[0] * 2 * M_PI;
  double k2 = xx[1] * 2 * M_PI;
  
  double kx = 0.5 * (k2 + k1);
  double ky = 0.5 * (k2 - k1);
  
  /* Sum over kz */
  for(int z=0; z < 2; z++){       
    double kz = M_PI * z;	  
    
    cx_double ek1 = ts()->ek1(kx, ky, kz);
    ff[0] += zk(ek1, ts()->tz, kz, delta());
  }
  
  return 0;   
}

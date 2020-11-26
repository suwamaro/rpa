/*****************************************************************************
*
* Functions for calculating the critical U.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "find_critical_U.h"

/* Member functions of FindUcIntegrandBilayer */
void FindUcIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double filling){
  assert(filling == 0.5);  
  hb_ = h;  
}

double FindUcIntegrandBilayer::integrand(const double *qvec) const {
  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  cx_double ek1 = ts()->ek1(kx, ky, kz);
  double n_minus = 1.0;  // Assume at a half filling
  double n_plus = 0.0;
  return (n_minus - n_plus) / std::abs(bk(up_spin, ek1, ts()->tz, kz));  // Spin does not matter.
}

int FindUcIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
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

    double qvec[3] = {kx, ky, kz};
    ff[0] += integrand(qvec);
  }
  
  return 0;   
}

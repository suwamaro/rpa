/*****************************************************************************
*
* Functions for finding the critical point.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "calc_single_particle_energy.h"
#include "find_critical_point.h"
#include "BinarySearch.h"

/* Member functions of FindQCPIntegrandBilayer */
void FindQCPIntegrand::set_mu(double _mu) { mu_ = _mu; }
hoppings2 *FindQCPIntegrand::ts() const { return ts_; }
double FindQCPIntegrand::mu() const { return mu_; }

/* Member functions of FindQCPIntegrandBilayer */
void FindQCPIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double _mu){
  hb_ = h;
  set_mu(_mu);
}

double FindQCPIntegrandBilayer::integrand(const double *qvec) const {
  /* Assume half filling */
  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  cx_double ek1 = ts()->ek1(kx, ky, kz);
  return 1. / std::abs(bk(up_spin, ek1, ts()->tz, kz));  // Spin does not matter.
}

int FindQCPIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
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

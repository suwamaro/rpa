/*****************************************************************************
*
* Functions for calculating the critical T.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "calc_single_particle_energy.h"
#include "find_critical_T.h"
#include "BinarySearch.h"

/* Member functions of FindTcIntegrandBilayer */
void FindTcIntegrand::set_T(double _T) { T_ = _T; }
void FindTcIntegrand::set_mu(double _mu) { mu_ = _mu; }
hoppings2 *FindTcIntegrand::ts() const { return ts_; }
double FindTcIntegrand::T() const { return T_; }
double FindTcIntegrand::mu() const { return mu_; }

/* Member functions of FindTcIntegrandBilayer */
void FindTcIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double _T, double _mu){
  hb_ = h;
  set_T(_T);
  set_mu(_mu);
}

double FindTcIntegrandBilayer::integrand(const double *qvec) const {
  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  cx_double ek1 = ts()->ek1(kx, ky, kz);

  /* Fermi density */
  double Em, Ep;
  std::tie(Em, Ep) = calc_single_particle_energy2(*ts(), kx, ky, kz, 0);  // delta == 0
  double n_minus = fermi_density(Em, kB*T(), mu());
  double n_plus = fermi_density(Ep, kB*T(), mu());	    	  
  return (n_minus - n_plus) / std::abs(bk(up_spin, ek1, ts()->tz, kz));  // Spin does not matter.
}

int FindTcIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
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

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
#include "calc_single_particle_energy.h"

/* Member functions of SelfConsistentIntegrand */
void SelfConsistentIntegrand::set_temperature(double _T){
  T_ = _T;
}

void SelfConsistentIntegrand::set_chemical_potential(double _mu){
  mu_ = _mu;
}

hoppings2 *SelfConsistentIntegrand::ts() const { return ts_; }
double SelfConsistentIntegrand::T() const { return T_; }
double SelfConsistentIntegrand::mu() const { return mu_; }


/* Member functions of SelfConsistentIntegrandBilayer */
void SelfConsistentIntegrandBilayer::set_parameters(hoppings_bilayer2 const& hb, double _mu, double _delta, double _T){
  set_temperature(_T);
  set_chemical_potential(_mu);
  hb_ = hb;
  delta_ = _delta;
}

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


/* Member functions of SelfConsistentIntegrand2 */
SelfConsistentIntegrand2::SelfConsistentIntegrand2(){
  max_iter_ = 200;
  eps_ = 1e-12;
  eps_func_ = 1e-10;
  non_zero_delta_lower_bound_ = 1e-12;  
}
std::size_t SelfConsistentIntegrand2::max_iter() const { return max_iter_; }
double SelfConsistentIntegrand2::eps() const { return eps_; }
double SelfConsistentIntegrand2::eps_func() const { return eps_func_; }
int SelfConsistentIntegrand2::L() const { return L_; }
double SelfConsistentIntegrand2::U() const { return U_; }
bool SelfConsistentIntegrand2::half_filling() const { return half_filling_; }
double SelfConsistentIntegrand2::filling() const { return filling_; }
bool SelfConsistentIntegrand2::T_equal_to_0() const { return T_equal_to_0_; }
double SelfConsistentIntegrand2::T() const { return T_; }
double SelfConsistentIntegrand2::delta() const { return delta_; }
double SelfConsistentIntegrand2::mu() const { return mu_; }
bool SelfConsistentIntegrand2::continuous_k() const { return continuous_k_; }
bool SelfConsistentIntegrand2::non_zero_delta() const { return non_zero_delta_; }
double SelfConsistentIntegrand2::non_zero_delta_lower_bound() const { return non_zero_delta_lower_bound_; }
double SelfConsistentIntegrand2::delta_upper_bound() const { return 0.5 * U(); }

void SelfConsistentIntegrand2::set_parameters(int _L, double _U, double _filling, double _T, double _delta, double _mu, bool _continuous_k, bool _non_zero_delta){
  L_ = _L;
  U_ = _U;
  filling_ = _filling;
  if (std::abs(filling() - 0.5) < 1e-12) {
    half_filling_ = true;    
  } else {
    half_filling_ = false;
  }
  T_ = _T;
  if (T() < 1e-15) {
    T_equal_to_0_ = true;
  } else {
    T_equal_to_0_ = false;
  }
  delta_ = _delta;
  mu_ = _mu;  
  continuous_k_ = _continuous_k;
  non_zero_delta_ = _non_zero_delta;
}

void SelfConsistentIntegrand2::set_input(double _delta, double _mu){
  delta_ = _delta;
  mu_ = _mu;
}

void SelfConsistentIntegrand2::set_max_iter(std::size_t _max_iter){
  max_iter_ = _max_iter;
}

void SelfConsistentIntegrand2::set_eps_func(double _eps){
  eps_func_ = _eps;
}

void SelfConsistentIntegrand2::set_mu_bounds(double lower, double upper){
  mu_lower_bound_ = lower;
  mu_upper_bound_ = upper;
}

bool SelfConsistentIntegrand2::invalid_params(double _delta, double _mu) const {
  if ( _delta < non_zero_delta_lower_bound() || _delta > delta_upper_bound() || _mu < mu_lower_bound() || _mu > mu_upper_bound() ) { return true; }
  else { return false; }
}

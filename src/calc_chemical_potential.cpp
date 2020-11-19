/*****************************************************************************
*
* Functions for calculating the chemical potential.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_chemical_potential.h"
#include "rpa_util.h"
#include "calc_single_particle_energy.h"

/* Member functions of ElecFillingIntegrand */
void ElecFillingIntegrand::set_parameters(double _mu, double _T, double _delta){
  mu_ = _mu;
  T_ = _T;
  delta_ = _delta;
}

/* Member functions of ElecFillingIntegrandBilayer */
void ElecFillingIntegrandBilayer::set_parameters(hoppings_bilayer2 const& hb, double _mu, double _T, double _delta){
  ElecFillingIntegrand::set_parameters(_mu, _T, _delta);
  hb_ = hb;
}

int ElecFillingIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
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

    /* Fermi density */
    double n_minus = 1.0;
    double n_plus = 0.0;
    if ( kB * T() < 1e-15 ) {
      n_minus = 1.0;
      n_plus = 0.0;
    } else {
      double Em, Ep;
      std::tie(Em, Ep) = calc_single_particle_energy2(*ts(), kx, ky, kz, delta());    
      n_minus = fermi_density(Em, kB*T(), mu());
      n_plus = fermi_density(Ep, kB*T(), mu());	    
    }
    
    ff[0] += n_minus + n_plus;
  }
  
  return 0;   
}


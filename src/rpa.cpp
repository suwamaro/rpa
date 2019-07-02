/*****************************************************************************
*
* RPA calculation for the fermionic Hubbard model
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"
#include "rpa_util.h"
#include "plot_chi0.h"
#include "self_consistent_eq.h"
#include "plot_self_consistent_eq_square.h"
#include "calc_gap.h"
#include "calc_dispersion.h"
#include "calc_velocity.h"
#include "plot_chi0_AF.h"
#include "calc_spectrum.h"

int main(){
  int L = 16;
  double eta = 0.00001;
    
  // plot_chi0( L, eta );
  // calc_size_dependence();
  // calc_chi();

  // double U = 5.0;  
  // double U = 0.337;
  double U = 0.34;
    
  // plot_self_consistent_eq_square( U );
  // plot_chi0_AF( U );
  // calc_dispersion_cubic( U, L );
  // calc_dispersion_square( U, L );
  
  // calc_spectrum_square( U, L, eta );
  // calc_spectrum_cubic( U, L, eta );
  
  // double theta = asin(1./sqrt(3));   // Cubic
  // double phi = 0;   // Octahedral rotation
  // double t3 = 0;   // Third-nearest xy-xy hopping amplitude
  
  // double theta = 0.234 * M_PI;   // Octahedral distortion
  double theta = 0.237 * M_PI;   // Octahedral distortion  
  double phi = 12. / 180. * M_PI;   // Octahedral rotation  
  double t3 = - 0.03;   // Third-nearest xy-xy hopping amplitude
  
  calc_spectrum_bilayer( theta, phi, t3, U, L, eta );    

  // calc_velocity_cubic();
  // calc_velocity_square();
  
  return 0;
}

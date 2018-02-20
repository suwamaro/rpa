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

int main(){
  // plot_chi0();
  // calc_size_dependence();
  // calc_chi();

  double U = 4.0;
  // plot_self_consistent_eq_square( U );

  // plot_chi0_AF( U );
  // calc_dispersion_cubic( U );
  calc_dispersion_square( U );

  // calc_velocity();
  
  return 0;
}

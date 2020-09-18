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
#include "find_critical_U.h"
#include "calc_spectrum.h"

int main(){
  /* Input parameters */
  rpa::parameters p("config.toml");

  /* Finding the critical U */
  find_critical_U_bilayer(p);
  
  /* Calculating the spectrum */
  calc_spectrum_bilayer2(p);
  
  return 0;
}

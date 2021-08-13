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
#include "find_critical_T.h"
#include "find_critical_point.h"
#include "calc_spectrum.h"
#include "calc_wave_func.h"
#include "calc_binding_energy.h"
#include "calc_phase_boundary.h"
#include "calc_current.h"

int main(int argc, char **argv){
  path base_dir;
  rpa::parameters p;
  std::tie(base_dir, p) = rpa::extract_parameters(argv[1]);
  
  /* Finding the critical U */
  find_critical_U_bilayer_output(base_dir, p);
  
  /* Finding the critical parameter */
  find_critical_point_bilayer_output(base_dir, p);

  /* Finding the critical T */
  find_critical_T_bilayer_output(base_dir, p);  

  /* Temperature dependence */
  solve_self_consistent_eqs_bilayer_T(base_dir, p);
  
  /* Calculating the spectrum */
  calc_spectrum_bilayer2(base_dir, p);

  /* Calculating the exciton wave function */
  calc_wave_func_bilayer(base_dir, p);

  /* Calculating the exciton binding energy */
  calc_binding_energy_bilayer(base_dir, p);  

  /* Obtaining the U-t4 (U-tz) phase diagram */
  calc_phase_boundary_U_bilayer(base_dir, p);
  calc_phase_boundary_t4_bilayer(base_dir, p);
  
  return 0;
}

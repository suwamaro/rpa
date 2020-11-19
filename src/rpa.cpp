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
#include "calc_spectrum.h"
#include "calc_wave_func.h"
#include "calc_binding_energy.h"

int main(int argc, char **argv){
  /* Simulation directory */
  path base_dir(argv[1]);
  
  /* Input parameters */  
  auto input_file = base_dir / "config.toml";
  if ( !exists(input_file) ) {
    std::cerr << "Required input file " << input_file << " does not exist.\n";
    std::exit(EXIT_FAILURE);
  }  
  rpa::parameters p(input_file.string());

  /* Finding the critical U */
  find_critical_U_bilayer(base_dir, p);

  /* Finding the critical U */
  find_critical_T_bilayer(base_dir, p);  
  
  /* Calculating the spectrum */
  calc_spectrum_bilayer2(base_dir, p);

  /* Calculating the exciton wave function */
  // calc_wave_func_bilayer(base_dir, p);

  /* Calculating the exciton binding energy */
  // calc_binding_energy_bilayer(base_dir, p);  
  
  return 0;
}

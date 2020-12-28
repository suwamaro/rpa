/*****************************************************************************
*
* Functions for calculating the critical U.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"
#include "parameters.h"

/* For bilayer lattices */
void calc_phase_boundary_U_bilayer(path& base_dir, rpa::parameters& pr);
void calc_phase_boundary_t4_bilayer(path& base_dir, rpa::parameters& pr);

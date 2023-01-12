/*****************************************************************************
*
* Functions for calculating the Raman scattering cross section.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _CALC_RAMAN_H
#define _CALC_RAMAN_H

#include "rpa.h"

void calc_Raman_bilayer(path& base_dir, rpa::parameters const& pr);
void calc_Raman_bilayer2(path& base_dir, rpa::parameters const& pr);

#endif // _CALC_RAMAN_H

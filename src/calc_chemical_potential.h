/*****************************************************************************
*
* Chemical potential calculation for the fermionic Hubbard model.
*
* Copyright (C) 2019 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __CALC_CHEMICAL_POTENTIAL__
#define __CALC_CHEMICAL_POTENTIAL__

#include "rpa_util.h"

double calc_chemical_potential_bilayer(int L, hoppings const& ts, double delta);
double calc_chemical_potential_bilayer2(int L, hoppings2 const& ts, double delta);

#endif // __CALC_CHEMICAL_POTENTIAL__

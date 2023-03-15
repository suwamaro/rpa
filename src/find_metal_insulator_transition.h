/*****************************************************************************
*
* Functions for finding a metal-insulator transition point.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "hoppings.h"

/* For bilayer lattices */
double find_metal_insulator_transition_t4_bilayer(rpa::parameters& pr);

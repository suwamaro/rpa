/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"

double self_consistent_eq_square(int L, double t, double mu, double delta);
double solve_self_consistent_eq_square(int L, double t, double mu, double U);

/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"

/* For a square lattice */
double self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double delta);
double solve_self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double U);

/* For a bilayer lattice */
double self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double mu, double delta);
double solve_self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double mu, double U);

/* For a simple cubic lattice */
double self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double delta);
double solve_self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double U);

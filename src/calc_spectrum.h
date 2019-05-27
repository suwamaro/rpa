/*****************************************************************************
*
* Functions for calculating the dispersion
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"

/* For a square lattice */
void calc_spectrum_square(double U, int L, double eta);

/* For a square lattice */
void calc_spectrum_bilayer(double theta, double phi, double t3, double U, int L, double eta);

/* For a simple cubic lattice */
void calc_spectrum_cubic(double U, int L, double eta);

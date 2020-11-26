/*****************************************************************************
*
* Functions for calculating the dispersion
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "parameters.h"

/* For square lattices */
void calc_spectrum_square(double U, int L, double eta);

/* For bilayer lattices */
void calc_spectrum_bilayer(double theta, double phi, double t3, double U, int L, double eta);
void calc_spectrum_bilayer2(path& base_dir, rpa::parameters const& pr);
extern "C" void calc_spectrum_bilayer2_wrapper(int argc, char **argv);

/* For simple cubic lattices */
void calc_spectrum_cubic(double U, int L, double eta);

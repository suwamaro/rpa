/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa.h"

/* Common to the lattices */
cx_double larger_eigenvalue(cx_double A, cx_double B, cx_double D);
double wave_vector_in_BZ(double k);
void add_to_sus_mat(cx_double& A, cx_double& B, cx_double& D, double e_free, double e_free2, double delta, cx_double omega);

/* For a square lattice */
double calc_eigval_square(int L, double t, double mu, double U, double delta, double qx, double qy, double omega);
double calc_gap_square(int L, double t, double mu, double U, double delta, double qx, double qy);
cx_double calc_intensity_square(int L, double t, double mu, double U, double delta, double qx, double qy, cx_double omega, int index);

/* For a simple cubic lattice */
double calc_eigval_cubic(int L, double t, double mu, double U, double delta, double qx, double qy, double qz, double omega);
double calc_gap_cubic(int L, double t, double mu, double U, double delta, double qx, double qy, double qz);

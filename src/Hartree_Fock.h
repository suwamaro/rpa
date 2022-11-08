/*****************************************************************************
*
* Functions for the Hartree-Fock approximation.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _HARTREE_FOCK_H
#define _HARTREE_FOCK_H

#include <armadillo>
#include "rpa.h"
#include "hoppings.h"
#include "parameters.h"
#include "cuba_helper.h"

double wave_vector_in_BZ(double k);
double BZ_factor(double kx, double ky);
int sublattice_spin_index(int g, int sigma);
int sign_spin_index(int g, int sigma);
int band_spin_index(int b, int sigma);
int spin_index(int sigma);
int index_to_spin(int idx);
int index_to_sign(int idx);
cx_vec gs_HF1(int spin, int sign, cx_double ek1, cx_double tz, double kz, double delta);
cx_mat gs_HF(cx_double ek1, cx_double tz, double kz, double delta);

# endif // _HARTREE_FOCK_H

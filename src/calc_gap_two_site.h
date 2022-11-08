/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _CALC_GAP_TWO_SITE_H
#define _CALC_GAP_TWO_SITE_H

#include <armadillo>
#include "rpa.h"
#include "hoppings.h"
#include "parameters.h"
#include "cuba_helper.h"
#include "Hartree_Fock.h"

std::tuple<double, double, double> calc_gap_two_site(hoppings_two_site const& ts, double mu, double U, double T, double delta, Polarization const& Pz, bool return_upper, bool verbose);

# endif // _CALC_GAP_TWO_SITE_H

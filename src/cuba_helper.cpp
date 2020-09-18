/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "cuba_helper.h"

/* Member functions of CubaParam */
CubaParam::CubaParam(rpa::parameters const& pr){
  epsrel = pr.epsrel;
  epsabs = pr.epsabs;
  flags = pr.flags;
  maxeval = pr.maxeval;
  key = pr.key;
}

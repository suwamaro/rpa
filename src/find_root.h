/*****************************************************************************
*
* Functions for finding a root.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __FIND_ROOT_H__
#define __FIND_ROOT_H__

#include <cmath>
#include <functional>
#include "rpa.h"

class NewtonRaphson {
public:
  explicit NewtonRaphson(int N);
  explicit NewtonRaphson(int N, int64_t nid);
  explicit NewtonRaphson(int N, int64_t nid, double mod_prefactor);
  virtual ~NewtonRaphson(){};
  double mod_factor(int64_t niter) const;
  vec calc_dx(int64_t niter) const;
  bool calc_dx2(int64_t niter, vec& dx) const;
  
  int64_t niter_decreasing;
  double mod_prefactor;
  vec F;  // Find a solution x such that F(x) == 0.
  mat J;  // Derivatives  
};

#endif  // __FIND_ROOT_H__

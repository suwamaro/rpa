/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __CUBA_HELPER__
#define __CUBA_HELPER__

#include "rpa.h"

/* Parameters for Cuba */
struct CubaParam {
  CubaParam(){};
  explicit CubaParam(rpa::parameters const& pr);
  ~CubaParam(){};

  const int NOPE = 2;  // +- and zz    
  const int NDIM = 2;  // 2D integral    
  const int NCOMP = NOPE * NSUBL*NSUBL * 2;  // 2 (+- and zz) * (NSUBL*NSUBL) (element) * 2 (real and imag)
  void *userdata = NULL;
  const int nvec = 1;
  double epsrel = 1e-5;
  double epsabs = 1e-10;
  int flags = 0;
  // int flags = 6;  // for debug
  int mineval = 0;
  int maxeval = 1 << 25;
  int key = 0; // cubature rule  
  const char *statefile = NULL;    
  void *spin = NULL;
};

# endif // __CUBA_HELPER__

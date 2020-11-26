#include <iostream>
#include "find_root.h"

/* Member functions of NewtonRaphson */
NewtonRaphson::NewtonRaphson(int N):F(N),J(N,N){
  niter_decreasing = 100;
}

NewtonRaphson::NewtonRaphson(int N, int64_t nid):niter_decreasing(nid),F(N),J(N,N){}

double NewtonRaphson::mod_factor(int64_t niter) const {
  /* 1 (niter < niter_decreasing) or ~1/niter (niter > niter_decreasing) */
  return niter > niter_decreasing ? 1. / (niter - niter_decreasing) : 1.0;
}

vec NewtonRaphson::calc_dx(int64_t niter){
  return mod_factor(niter) * arma::solve(J, - F);
}

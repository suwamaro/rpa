#include <iostream>
#include "find_root.h"

/* Member functions of NewtonRaphson */
NewtonRaphson::NewtonRaphson(int N):F(N),J(N,N){
  niter_decreasing = 100;
  mod_prefactor = 1.0;
}

NewtonRaphson::NewtonRaphson(int N, int64_t nid):niter_decreasing(nid),F(N),J(N,N){
  mod_prefactor = 1.0;
}

NewtonRaphson::NewtonRaphson(int N, int64_t nid, double mod_pf):niter_decreasing(nid),mod_prefactor(mod_pf),F(N),J(N,N){}

double NewtonRaphson::mod_factor(int64_t niter) const {
  /* 1 (niter < niter_decreasing) or ~1/niter (niter > niter_decreasing) */
  double factor = niter > niter_decreasing ? 1. / (niter - niter_decreasing) : 1.0;
  return mod_prefactor * factor;
}

vec NewtonRaphson::calc_dx(int64_t niter) const {
  return mod_factor(niter) * arma::solve(J, - F);
}

bool NewtonRaphson::calc_dx2(int64_t niter, vec& dx) const {
  bool status = arma::solve(dx, J, - F);
  dx *= mod_factor(niter);
  return status;
}

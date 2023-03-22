/*****************************************************************************
*
* Functions for finding a root.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <iostream>
#include "optimization.h"

/* Member functions of NewtonRaphson */
NelderMead::NelderMead(int d):dim(d),t(0),xs(d+1),fs(d+1),alpha(1.),gamma(2.),rho(0.5),sigma(0.5){
  init();
}
NelderMead::NelderMead(int d, double _alpha, double _gamma, double _rho, double _sigma):dim(d),t(0),xs(d+1),fs(d+1),alpha(_alpha),gamma(_gamma),rho(_rho),sigma(_sigma){
  init();
}

void NelderMead::init(){
  x0.set_size(dim);
  eps = 1e-10;
}

void NelderMead::init_x(std::vector<vec> const& _xs){
  xs = _xs;
}

void NelderMead::set_eps(double _eps){
  eps = _eps;
}

void NelderMead::reset(){
  t = 0;
}

void NelderMead::sort(){
  std::sort(xs.begin(), xs.end(), compare);
  for(int i=0; i <= dim; ++i){
    fs[i] = f(xs[i]);
  }
  is_sorted = true;
}

void NelderMead::set_x0(){
  x0.zeros();
  for(int i=0; i < dim; ++i){
    x0 += xs[i];
  }
  x0 /= (double)(dim);
}

void NelderMead::step(){
  ++t;
  
  /* Order */
  if (!is_sorted) { sort(); }

  /* Setting is_sorted false. */
  is_sorted = false;
  
  /* Calculate x0 */
  set_x0();
  
  /* Reflection */
  vec xr = x0 + alpha * (x0 - xs[dim]);
  double fxr = f(xr);
  
  if (fs[0] <= fxr && fxr <= fs[dim-1]) {
    xs[dim] = xr;
    return;
  }

  /* Expansion */
  if (fxr < fs[0]){
    vec xe = x0 + gamma * (xr - x0);
    double fxe = f(xe);
    if (fxe < fxr) {
      xs[dim] = xe;
    } else {
      xs[dim] = xr;
    }
    return;
  }

  /* Contraction */
  if (fxr < fs[dim]) {
    vec xc = x0 + rho * (xr - x0);
    double fxc = f(xc);
    if (fxc < fxr) {
      xs[dim] = xc;
      return;
    }
  } else {
    vec xc = x0 + rho * (xs[dim] - x0);
    double fxc = f(xc);
    if (fxc < fs[dim]) {
      xs[dim] = xc;
      return;
    }
  }

  /* Shrink */
  for(int i=1; i <= dim; ++i){
    xs[i] = xs[0] + sigma * (xs[i] - xs[0]);
  }

  return;
}

double NelderMead::distance(vec const& x1, vec const& x2) const {
  double sum = 0.;
  for(int i=0; i < dim; ++i){
    double diff = x2[i] - x1[i];
    sum += diff * diff;
  }
  return sum;
}

bool NelderMead::is_terminated(){
  set_x0();
  for(int i=0; i <= dim; ++i){
    double di = distance(xs[i], x0);
    if (di > eps * eps) {
      return false;
    }
  }
  return true;
}

void NelderMead::output(ostream& out) const {
  out << t << " ";
  for(vec x: xs){
    for(int d=0; d < dim; ++d){
      out << "  " << x[d];
    }
  }
  out << "  ";
  for(double fx: fs){
    out << "  " << fx;
  }
  out << std::endl;
}

void NelderMead::get_result(vec& x_opt, double& f_opt){
  sort();
  x_opt = xs[0];
  f_opt = fs[0];
}

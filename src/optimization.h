/*****************************************************************************
*
* Functions for optimization.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __OPTIMIZATION_H__
#define __OPTIMIZATION_H__

#include <cmath>
//#include <functional>
#include "rpa.h"

class NelderMead {
public:
  explicit NelderMead(int d);
  explicit NelderMead(int d, double alpha, double gamma, double rho, double sigma);  
  virtual ~NelderMead(){};

  void init();  
  void init_x(std::vector<vec> const& _xs);
  void reset();
  void sort();
  void set_x0();
  void step();
  double (*f)(vec const&);
  bool (*compare)(vec const& x1, vec const& x2);  
  bool is_terminated() const;
  void set_tolerance(double tol);
  void output(ostream& out) const;
  void get_result(vec& x_opt, double& f_opt);

private:
  int dim;
  std::size_t t;
  double tolerance = 1e-12;
  bool is_sorted;
  std::vector<vec> xs;
  std::vector<double> fs;
  vec x0;
  double alpha;
  double gamma;
  double rho;
  double sigma;
};

#endif  // __OPTIMIZATION_H__

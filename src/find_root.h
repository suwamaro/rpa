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

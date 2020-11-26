#include <cmath>
#include <functional>
#include "rpa.h"

class NewtonRaphson {
public:
  explicit NewtonRaphson(int N);
  explicit NewtonRaphson(int N, int64_t nid);
  virtual ~NewtonRaphson(){};
  double mod_factor(int64_t niter) const;
  vec calc_dx(int64_t niter);

  int64_t niter_decreasing;
  vec F;  // Find a solution x such that F(x) == 0.
  mat J;  // Derivatives  
};

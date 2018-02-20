#include <cmath>
#include <functional>

class BinarySearch {
private:
  unsigned int max_iter = 1e+8;
  double eps_fx = 1e-10;
  
public:
  BinarySearch(){};
  void find_solution(double& x, double target, std::function<double(double x)> const& f, bool additive = false, double x_delta = 0.01, double x_MIN = NAN, double x_MAX = NAN, bool debug = false);
  
  virtual ~BinarySearch(){};
};

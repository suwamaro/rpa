#include <cmath>
#include <functional>

class BinarySearch {
public:
  BinarySearch(){};
  virtual ~BinarySearch(){};
  void set_x_MIN(double _x_MIN) { x_MIN_ = _x_MIN; }
  void set_x_MAX(double _x_MAX) { x_MAX_ = _x_MAX; }  
  bool find_solution(double& x, double target, std::function<double(double x)> const& f, bool additive = false, double x_delta = 0.01, bool debug = false);

private:
  unsigned int max_iter = 50;
  double eps_fx = 1e-10;
  double x_MIN_ = NAN;
  double x_MAX_ = NAN;  
};

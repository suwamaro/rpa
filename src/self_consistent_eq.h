/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"
#include "find_root.h"

/* Number of parameters for solving the self-consistent equations. */
constexpr int n_sc_params = 2;

/* Integrand */
class SelfConsistentIntegrand {
public:
  explicit SelfConsistentIntegrand(hoppings2 *ts):ts_(ts){}
  virtual ~SelfConsistentIntegrand(){}
  virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
  void set_temperature(double _T);
  void set_chemical_potential(double _mu);
  hoppings2 *ts() const;
  double T() const;
  double mu() const;
  
private:  
  hoppings2 *ts_;
  double T_;
  double mu_;
};

class SelfConsistentIntegrandBilayer : public SelfConsistentIntegrand {
public:
  SelfConsistentIntegrandBilayer():SelfConsistentIntegrand(&hb_){}
  virtual ~SelfConsistentIntegrandBilayer(){}
  void set_parameters(hoppings_bilayer2 const& ts, double _mu, double _delta, double _T);
  void set_delta(double _delta);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
  double delta() const { return delta_; }
  
private:
  hoppings_bilayer2 hb_;
  double delta_;  
};

class SelfConsistentIntegrand2 {
public:
  explicit SelfConsistentIntegrand2(hoppings2 *ts);
  virtual ~SelfConsistentIntegrand2(){}
  virtual void set_parameters(int _L, double _U, double _filling, double _T, double _delta, double _mu, bool _continuous_k, bool _non_zero_delta);
  void set_input(double _delta, double _mu);
  void set_eps_func(double _eps);
  int64_t max_iter() const;
  double eps() const;
  double eps_func() const;
  hoppings2 *ts() const;
  int L() const;
  double U() const;
  double filling() const;
  double T() const;
  double delta() const;
  double mu() const;
  bool continuous_k() const;
  bool non_zero_delta() const;
  double non_zero_delta_lower_bound() const;
  double delta_upper_bound() const;
  
private:
  int64_t max_iter_;
  double eps_;
  double eps_func_;
  hoppings2 *ts_;
  int L_;
  double U_;
  double filling_;
  double T_;
  double delta_;
  double mu_;
  bool continuous_k_;
  bool non_zero_delta_;
  double non_zero_delta_lower_bound_;
};

class SelfConsistentIntegrand2Bilayer : public SelfConsistentIntegrand2 {
public:
  SelfConsistentIntegrand2Bilayer();
  virtual ~SelfConsistentIntegrand2Bilayer(){}
  void set_parameters(rpa::parameters const& pr, int _L, hoppings_bilayer2 const& ts, double _U, double _filling, double _T, double _delta, double _mu, bool _continous_k, bool _non_zero_delta);
  std::tuple<double, double> calc_elec_density(const double *qvec) const;
  std::tuple<double, double> calc_compressibility(const double *qvec) const;
  double integrand_mean_field(const double *qvec) const;
  double integrand_elec_density(const double *qvec) const;
  double integrand_mean_field_der_delta(const double *qvec) const;
  double integrand_mean_field_der_mu(const double *qvec) const;
  double integrand_elec_density_der_delta(const double *qvec) const;
  double integrand_elec_density_der_mu(const double *qvec) const;
  std::tuple<double, double> calc_qs(const cubareal *xx) const;
  void integral(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
  double calc() const;
  double calc_mean_field();
  double calc_elec_density();
  double calc_mean_field_der_delta();
  double calc_mean_field_der_mu();
  double calc_elec_density_der_delta();
  double calc_elec_density_der_mu();
  void update_parameters(int64_t niter, double& delta, double& mu);  
  bool find_solution(double& delta, double& mu, bool verbose);
  
private:
  hoppings_bilayer2 hb_;
  CubaParam cbp_;
  NewtonRaphson nr_;
  double (SelfConsistentIntegrand2Bilayer::*integrand_ptr)(const double*) const;
};


/* For a square lattice */
double self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double delta);
double solve_self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double U);

/* For a bilayer lattice */
double self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double delta);
double solve_self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double U);
double solve_self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double U, double filling, double T, CubaParam const& cbp, bool continuous_k);
/* For finite temperatures */
std::tuple<double, double> solve_self_consistent_eqs_bilayer(rpa::parameters const& pr, int L, hoppings_bilayer2 const& ts, double U, double filling, double T, bool continuous_k, bool non_zero_delta);
void solve_self_consistent_eqs_bilayer_T(path& base_dir, rpa::parameters const& pr);

/* For a simple cubic lattice */
double self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double delta);
double solve_self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double U);

/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <gsl/gsl_multiroots.h>
#include <cuba.h>
#include "rpa.h"
#include "hoppings.h"
#include "cuba_helper.h"
#include "find_root.h"
#include "optimization.h"

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
  SelfConsistentIntegrand2();
  virtual ~SelfConsistentIntegrand2(){}
  virtual void set_parameters(int _L, double _U, double _filling, double _T, double _delta, double _mu, bool _continuous_k, bool _non_zero_delta);
  void set_max_iter(std::size_t max_iter);
  void set_input(double _delta, double _mu);
  void set_U(double U);
  void set_eps_func(double _eps);
  std::size_t max_iter() const;
  int precision() const { return precision_; }
  int precision2() const { return precision2_; }  
  double eps() const;
  double eps_func() const;
  int L() const;
  double U() const;
  bool half_filling() const;
  double filling() const;
  bool T_equal_to_0() const;
  double T() const;
  double delta() const;
  double mu() const;
  bool continuous_k() const;
  bool non_zero_delta() const;
  double non_zero_delta_lower_bound() const;
  double delta_upper_bound() const;
  void set_mu_bounds(double lower, double upper);
  double mu_lower_bound() const { return mu_lower_bound_; }
  double mu_upper_bound() const { return mu_upper_bound_; }
  virtual bool invalid_params(double _delta, double _mu) const;
  
private:
  std::size_t max_iter_;
  int precision_ = 12;
  int precision2_ = 15;  
  double eps_;
  double eps_func_;
  int L_;
  double U_;
  bool half_filling_;
  double filling_;
  bool T_equal_to_0_;
  double T_;
  double delta_;
  double mu_;
  bool continuous_k_;
  bool non_zero_delta_;
  double non_zero_delta_lower_bound_;
  double mu_lower_bound_ = std::numeric_limits<double>::min();
  double mu_upper_bound_ = std::numeric_limits<double>::max();
};

class SelfConsistentIntegrand2Bilayer : public SelfConsistentIntegrand2 {
public:
  SelfConsistentIntegrand2Bilayer();
  virtual ~SelfConsistentIntegrand2Bilayer(){}
  void set_parameters(rpa::parameters const& pr, int _L, hoppings_bilayer2 const& ts, double _U, double _filling, double _T, double _delta, double _mu, bool _continous_k, bool _non_zero_delta);
  const hoppings_bilayer2 *ts() const;
  std::tuple<double, double> calc_elec_density(const double *qvec) const;
  std::tuple<double, double> calc_compressibility(const double *qvec) const;
  double calc_energy(const double *qvec) const;    
  double integrand_mean_field(const double *qvec) const;
  double integrand_elec_density(const double *qvec) const;
  double integrand_mean_field_der_delta(const double *qvec) const;
  double integrand_mean_field_der_mu(const double *qvec) const;
  double integrand_elec_density_der_delta(const double *qvec) const;
  double integrand_elec_density_der_mu(const double *qvec) const;
  double integrand_energy(const double *qvec) const;  
  std::tuple<double, double> calc_qs(const cubareal *xx) const;
  void integral(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
  double calc() const;
  double calc_mean_field_function();  
  double calc_mean_field();
  double calc_elec_density();
  double calc_mean_field_der_delta();
  double calc_mean_field_der_mu();
  double calc_elec_density_der_delta();
  double calc_elec_density_der_mu();
  double calc_energy();
  double calc_diff_nr();
  double calc_diff();  
  double calc_diff(vec const& x);
  double calc_diff(double delta, double mu);
  double calc_chemical_potential(double delta, bool verbose=false) const;  
  void update_parameters(std::size_t niter, double& delta, double& mu);
  int gsl_function_helper(const gsl_vector *x, void *params, gsl_vector *f);
  int gsl_function_df_helper(const gsl_vector *x, void *params, gsl_matrix *J);
  int gsl_function_fdf_helper(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);
  bool find_solution_gsl(double& delta, double& mu, bool verbose);
  bool find_solution_nr(double& delta, double& mu, bool verbose);
  bool find_solution_nm(double& delta, double& mu, bool verbose);
  bool find_solution_using_1d_solver(double& delta, double& mu, bool verbose);
  bool find_solution_default(double& delta, double& mu, bool verbose);    
  bool find_solution(double& delta, double& mu, bool verbose);
  double NelderMead_f(vec const& x) const { return nm_.f(x); }
  bool use_gsl() const { return use_gsl_; }
  bool use_gsl_fdf_solver() const { return use_gsl_fdf_solver_; }
  bool use_NewtonRaphson() const { return use_NewtonRaphson_; }  
  bool use_NelderMead() const { return use_NelderMead_; }
  bool use_1d_solver() const { return use_1d_solver_; }  
  
private:
  int n_sc_params = 2;  
  hoppings_bilayer2 hb_;
  CubaParam cbp_;
  NewtonRaphson nr_;
  NelderMead nm_;
  bool use_gsl_;
  bool use_gsl_fdf_solver_;
  bool use_NewtonRaphson_;  
  bool use_NelderMead_;
  bool use_1d_solver_;
  double (SelfConsistentIntegrand2Bilayer::*integrand_ptr)(const double*) const;
};


/* For a square lattice */
double self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double delta);
double solve_self_consistent_eq_square(int L, hoppings_square const& ts, double mu, double U);

/* For a bilayer lattice */
double self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double delta);
double solve_self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double U);
double solve_self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double U, double mu, double T, CubaParam const& cbp, bool continuous_k);
/* For finite temperatures */
std::tuple<double, double> solve_self_consistent_eqs_bilayer(rpa::parameters const& pr, int L, hoppings_bilayer2 const& ts, double U, double filling, double T, bool continuous_k, double delta_i = 0.01, double mu_i = 0.);
std::tuple<double, double> solve_self_consistent_eqs_bilayer2(rpa::parameters const& pr, int L, hoppings_bilayer2 const& ts, double U, double filling, double T, bool continuous_k, bool non_zero_delta, double delta_i = 0.01, double mu_i = 0.);   // with non_zero_delta
void solve_self_consistent_eqs_bilayer_T(path& base_dir, rpa::parameters const& pr);

/* For a simple cubic lattice */
double self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double delta);
double solve_self_consistent_eq_cubic(int L, hoppings_cubic const& ts, double mu, double U);

/* For the two-site problem */
double solve_self_consistent_eq_two_site(hoppings_two_site const& ts, double U, double T);

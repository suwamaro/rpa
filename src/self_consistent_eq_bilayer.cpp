/*****************************************************************************
*
* Functions for solving the self-consistent equation
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "self_consistent_eq.h"
#include "rpa_util.h"
#include "calc_gap.h"
#include "calc_single_particle_energy.h"
#include "calc_chemical_potential.h"
#include "BinarySearch.h"
#include "find_critical_U.h"
#include "find_critical_T.h"

/* Creating an instance */
SelfConsistentIntegrandBilayer scib;
SelfConsistentIntegrand2Bilayer sci2b;

int sci2b_gsl_function(const gsl_vector *x, void *params, gsl_vector *f){
  return sci2b.gsl_function_helper(x, params, f);
}
int sci2b_gsl_function_df(const gsl_vector *x, void *params, gsl_matrix *J){
  return sci2b.gsl_function_df_helper(x, params, J);
}
int sci2b_gsl_function_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J){
  return sci2b.gsl_function_fdf_helper(x, params, f, J);
}
  
double sci2b_calc_diff(vec const& x){
  return sci2b.calc_diff(x);
}

bool compare_x(vec const& x1, vec const& x2){
  double f1 = sci2b.NelderMead_f(x1);
  double f2 = sci2b.NelderMead_f(x2);
  return f1 < f2;
}

/* Wrapper functions */
int sc_integrand_bilayer(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return scib.calc(ndim, xx, ncomp, ff, userdata);
}

int sc_integrand_bilayer2(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  sci2b.integral(ndim, xx, ncomp, ff, userdata);
  return 0;
}

/* Member functions of SelfConsistentIntegrand2Bilayer */
SelfConsistentIntegrand2Bilayer::SelfConsistentIntegrand2Bilayer():SelfConsistentIntegrand2(),nr_(n_sc_params),nm_(n_sc_params),use_gsl_(false),use_gsl_fdf_solver_(false),use_NewtonRaphson_(false),use_NelderMead_(false),use_1d_solver_(false){
  nm_.f = &sci2b_calc_diff;
  nm_.compare = &compare_x;
}

void SelfConsistentIntegrand2Bilayer::set_parameters(rpa::parameters const& pr, int _L, hoppings_bilayer2 const& _hb, double _U, double _filling, double _T, double _delta, double _mu, bool _continuous_k, bool _non_zero_delta){
  SelfConsistentIntegrand2::set_parameters(_L, _U, _filling, _T, _delta, _mu, _continuous_k, _non_zero_delta);
  SelfConsistentIntegrand2::set_max_iter(pr.max_iter);  
  SelfConsistentIntegrand2::set_eps_func(pr.epsfunc);
  hb_ = _hb;
  cbp_.set_parameters(pr);
  nr_.mod_prefactor = pr.mod_prefactor;
  use_NewtonRaphson_ = pr.use_NewtonRaphson;  
  use_NelderMead_ = pr.use_NelderMead;
  use_1d_solver_ = pr.use_1d_solver;  
}

const hoppings_bilayer2 *SelfConsistentIntegrand2Bilayer::ts() const { return &hb_; }

std::tuple<double, double> SelfConsistentIntegrand2Bilayer::calc_elec_density(const double *qvec) const {
  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  double Em, Ep;
  std::tie(Em, Ep) = calc_single_particle_energy2(*ts(), kx, ky, kz, delta());
  double n_minus = fermi_density(Em, kB*T(), mu());
  double n_plus = fermi_density(Ep, kB*T(), mu());
  return std::make_tuple(n_minus, n_plus);
}

std::tuple<double, double> SelfConsistentIntegrand2Bilayer::calc_compressibility(const double *qvec) const {
  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  double Em, Ep;
  std::tie(Em, Ep) = calc_single_particle_energy2(*ts(), kx, ky, kz, delta());
  double g_minus = compressibility(Em, kB*T(), mu());
  double g_plus = compressibility(Ep, kB*T(), mu());
  return std::make_tuple(g_minus, g_plus);
}

double SelfConsistentIntegrand2Bilayer::calc_energy(const double *qvec) const {
  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];  
  double Em, Ep;
  std::tie(Em, Ep) = calc_single_particle_energy2(*ts(), kx, ky, kz, delta());
  double f_minus = fermi_energy(Em, kB*T(), mu());
  double f_plus = fermi_energy(Ep, kB*T(), mu());
  return f_minus + f_plus;   // Free energy
}

double SelfConsistentIntegrand2Bilayer::integrand_mean_field(const double *qvec) const {
  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  cx_double ek1 = ts()->ek1(kx, ky, kz);
  double n_minus, n_plus;
  std::tie(n_minus, n_plus) = calc_elec_density(qvec);
  return (n_minus - n_plus) * zk_over_delta(ek1, ts()->tz, kz, delta());  // Using zk_over_delta
}

double SelfConsistentIntegrand2Bilayer::integrand_elec_density(const double *qvec) const {
  double n_minus, n_plus;
  std::tie(n_minus, n_plus) = calc_elec_density(qvec);
  return n_minus + n_plus;
}

double SelfConsistentIntegrand2Bilayer::integrand_mean_field_der_delta(const double *qvec) const {
  double n_minus, n_plus;
  std::tie(n_minus, n_plus) = calc_elec_density(qvec);
  
  double g_minus, g_plus;
  std::tie(g_minus, g_plus) = calc_compressibility(qvec);

  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  cx_double ek1 = ts()->ek1(kx, ky, kz);
  double zod = zk_over_delta(ek1, ts()->tz, kz, delta());
  return ((g_minus + g_plus) * std::pow(zod, 2) - (n_minus - n_plus) * std::pow(zod, 3)) * delta();
}

double SelfConsistentIntegrand2Bilayer::integrand_mean_field_der_mu(const double *qvec) const {
  double g_minus, g_plus;
  std::tie(g_minus, g_plus) = calc_compressibility(qvec);

  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  cx_double ek1 = ts()->ek1(kx, ky, kz);
  double zod = zk_over_delta(ek1, ts()->tz, kz, delta());
  return (g_minus - g_plus) * zod;
}

double SelfConsistentIntegrand2Bilayer::integrand_elec_density_der_delta(const double *qvec) const {
  double g_minus, g_plus;
  std::tie(g_minus, g_plus) = calc_compressibility(qvec);

  double kx = qvec[0];
  double ky = qvec[1];
  double kz = qvec[2];
  cx_double ek1 = ts()->ek1(kx, ky, kz);
  double z = zk(ek1, ts()->tz, kz, delta());
  return (g_minus - g_plus) * z;  
}

double SelfConsistentIntegrand2Bilayer::integrand_elec_density_der_mu(const double *qvec) const {
  double g_minus, g_plus;
  std::tie(g_minus, g_plus) = calc_compressibility(qvec);
  return g_minus + g_plus;
}

double SelfConsistentIntegrand2Bilayer::integrand_energy(const double *qvec) const {
  return calc_energy(qvec);
}

std::tuple<double, double> SelfConsistentIntegrand2Bilayer::calc_qs(const cubareal *xx) const {
  /* Wavenumbers */
  double k1 = xx[0] * 2 * M_PI;
  double k2 = xx[1] * 2 * M_PI;
  
  double kx = 0.5 * (k2 + k1);
  double ky = 0.5 * (k2 - k1);
  return std::make_tuple(kx, ky);
}

void SelfConsistentIntegrand2Bilayer::integral(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {  
  /* Reset */
  ff[0] = 0;
  
  double kx, ky;
  std::tie(kx, ky) = calc_qs(xx);
  
  /* Sum over kz */
  for(int z=0; z < 2; ++z){       
    double kz = M_PI * z;
    double qvec[] = {kx, ky, kz};
    ff[0] += (this->*integrand_ptr)(qvec);    
  }
}

double SelfConsistentIntegrand2Bilayer::calc_mean_field_function(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_mean_field;
  return calc() * delta();
}

double SelfConsistentIntegrand2Bilayer::calc_mean_field(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_mean_field;
  return calc() - 1. / U();  // == 0 for solutions
}

double SelfConsistentIntegrand2Bilayer::calc_elec_density(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_elec_density;
  return calc() - filling();  // == 0 for solutions
}

double SelfConsistentIntegrand2Bilayer::calc_mean_field_der_delta(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_mean_field_der_delta;
  return calc();
}

double SelfConsistentIntegrand2Bilayer::calc_mean_field_der_mu(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_mean_field_der_mu;
  return calc();
}

double SelfConsistentIntegrand2Bilayer::calc_elec_density_der_delta(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_elec_density_der_delta;
  return calc();
}

double SelfConsistentIntegrand2Bilayer::calc_elec_density_der_mu(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_elec_density_der_mu;
  return calc();
}

double SelfConsistentIntegrand2Bilayer::calc_energy(){
  integrand_ptr = &SelfConsistentIntegrand2Bilayer::integrand_energy;
  return calc();
}

double SelfConsistentIntegrand2Bilayer::calc_diff_nr(){
  if ( non_zero_delta() ) {
    return std::norm(nr_.F(0)) + std::norm(nr_.F(1));
  } else {
    return std::norm(nr_.F(1));
  }
}

double SelfConsistentIntegrand2Bilayer::calc_diff(){
  double f1 = calc_elec_density();
  if ( non_zero_delta() ) {
    double f0 = calc_mean_field();    
    return f0 * f0 + f1 * f1;
  } else {
    return f1 * f1;
  }
}

double SelfConsistentIntegrand2Bilayer::calc_diff(vec const& x){
  set_input(x(0), x(1));
  return calc_diff();
}

double SelfConsistentIntegrand2Bilayer::calc_diff(double _delta, double _mu){
  set_input(_delta, _mu);
  return calc_diff();
}

void SelfConsistentIntegrand2Bilayer::update_parameters(std::size_t niter, double& _delta, double& _mu){
  if ( non_zero_delta() ) {
    /* Updating delta and mu */
    nr_.J(0,0) = calc_mean_field_der_delta();
    nr_.J(0,1) = calc_mean_field_der_mu();
    nr_.J(1,0) = calc_elec_density_der_delta();
    nr_.J(1,1) = calc_elec_density_der_mu();
    nr_.F(0) = calc_mean_field();
    nr_.F(1) = calc_elec_density();
    vec dx;
    bool status = nr_.calc_dx2(niter, dx);
    _delta += dx(0);
    if ( _delta < non_zero_delta_lower_bound() ) { _delta = non_zero_delta_lower_bound(); }
    if ( _delta > delta_upper_bound() ) { _delta = delta_upper_bound(); }    
    _mu += dx(1);
  } else {
    /* Updating mu only */
    _delta = 0;
    nr_.J(1,1) = calc_elec_density_der_mu();
    nr_.F(1) = calc_elec_density();
    _mu += - nr_.mod_factor(niter) * nr_.F(1) / nr_.J(1,1);
  }
  if ( _mu < mu_lower_bound() ) { _mu = mu_lower_bound(); }
  if ( _mu > mu_upper_bound() ) { _mu = mu_upper_bound(); }
}

int SelfConsistentIntegrand2Bilayer::gsl_function_helper(const gsl_vector *x, void *params, gsl_vector *f){  
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  if (invalid_params(x0, x1)) {
    gsl_vector_set(f, 0, 1e+12);
    gsl_vector_set(f, 1, 1e+12);
  } else {    
    set_input(x0, x1);
    const double y0 = calc_mean_field();
    const double y1 = calc_elec_density();
    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
  }
  
  return GSL_SUCCESS;
}

int SelfConsistentIntegrand2Bilayer::gsl_function_df_helper(const gsl_vector *x, void *params, gsl_matrix *J){
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);
  
  set_input(x0, x1);  
  const double df00 = calc_mean_field_der_delta();
  const double df01 = calc_mean_field_der_mu();
  const double df10 = calc_elec_density_der_delta();
  const double df11 = calc_elec_density_der_mu();
  gsl_matrix_set(J, 0, 0, df00);
  gsl_matrix_set(J, 0, 1, df01);
  gsl_matrix_set(J, 1, 0, df10);
  gsl_matrix_set(J, 1, 1, df11);
  
  return GSL_SUCCESS;
}

int SelfConsistentIntegrand2Bilayer::gsl_function_fdf_helper(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J){
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  if (invalid_params(x0, x1)) {
    gsl_vector_set(f, 0, 1e+12);
    gsl_vector_set(f, 1, 1e+12);
  } else {      
    set_input(x0, x1);
  
    const double y0 = calc_mean_field();
    const double y1 = calc_elec_density();  
    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);

    const double df00 = calc_mean_field_der_delta();
    const double df01 = calc_mean_field_der_mu();
    const double df10 = calc_elec_density_der_delta();
    const double df11 = calc_elec_density_der_mu();
    gsl_matrix_set(J, 0, 0, df00);
    gsl_matrix_set(J, 0, 1, df01);
    gsl_matrix_set(J, 1, 0, df10);
    gsl_matrix_set(J, 1, 1, df11);
  }
  
  return GSL_SUCCESS;  
}

bool SelfConsistentIntegrand2Bilayer::find_solution_gsl(double& _delta, double& _mu, bool verbose){
  const gsl_multiroot_fsolver_type *T;
  const gsl_multiroot_fdfsolver_type *Tdf;        
  gsl_multiroot_fsolver *s;
  gsl_multiroot_fdfsolver *sdf;

  int status;
  std::size_t iter = 0;

  const size_t n = 2;
  gsl_multiroot_function f = {&sci2b_gsl_function, n, NULL};
  gsl_multiroot_function_fdf fdf = {&sci2b_gsl_function, &sci2b_gsl_function_df, &sci2b_gsl_function_fdf, n, NULL};    

  set_input(_delta, _mu);
    
  double x_init[2] = {delta(), mu()};
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, x_init[0]);
  gsl_vector_set(x, 1, x_init[1]);

  T = gsl_multiroot_fsolver_hybrids;    
  Tdf = gsl_multiroot_fdfsolver_hybridsj;
    
  s = gsl_multiroot_fsolver_alloc(T, n);
  sdf = gsl_multiroot_fdfsolver_alloc(Tdf, n);
    
  gsl_multiroot_fsolver_set(s, &f, x);    
  gsl_multiroot_fdfsolver_set(sdf, &fdf, x);

  /* Iteration */
  do {
    iter++;
    if (use_gsl_fdf_solver()) {
      status = gsl_multiroot_fdfsolver_iterate(sdf);
    } else {
      status = gsl_multiroot_fsolver_iterate(s);
    }
    if (status) break;
    status = gsl_multiroot_test_residual(s->f, eps_func());
  } while(status == GSL_CONTINUE && iter <= max_iter());

  printf ("status = %s\n", gsl_strerror(status));
  printf ("number of iterations = %zu\n", iter);
  printf ("root found at (%g, %g)\n",
	  gsl_vector_get(s->x, 0),
	  gsl_vector_get(s->x, 1));

  /* Getting the result. */
  _delta = gsl_vector_get(x, 0);
  _mu = gsl_vector_get(x, 1);

  /* Free memory */
  gsl_multiroot_fsolver_free(s);
  gsl_multiroot_fdfsolver_free(sdf);    
  gsl_vector_free(x);

  return true;  
}

bool SelfConsistentIntegrand2Bilayer::find_solution_nr(double& _delta, double& _mu, bool verbose){
  std::cout << "Using the Newton-Raphson method." << std::endl;
  double delta2 = _delta;
  double mu2 = _mu;
  std::size_t niter = 0;
  double diff = 0;
  do {
    ++niter;
    if ( niter == max_iter() ) {
      std::cout << "Number of iteration reaches the limit " << max_iter() << std::endl;
      std::cout << "The remained error squared is " << std::setprecision(precision()) << diff << std::endl;
      break;
    }
    _delta = delta2;
    _mu = mu2;
    set_input(_delta, _mu);    
    update_parameters(niter, delta2, mu2);
    diff = calc_diff_nr();
    if ( verbose ) {
      if ( non_zero_delta() ) {
	std::cout << niter << "  "
		  << _delta << "  " << delta2 << "  " << _mu << "  " << mu2 << "  "
		  << diff << "  "
		  << nr_.F(0) << "  " << nr_.F(1) << "  "
		  << nr_.J(0,0) << "  " << nr_.J(0,1) << "  " << nr_.J(1,0) << "  " << nr_.J(1,1) << std::endl;
      } else {
	std::cout << niter << "  " << _delta << "  "
		  << _mu << "  " << mu2 << "  "
		  << diff << "  "
		  << nr_.F(1) << "  "
		  << nr_.J(1,1)
		  << std::endl;
      }
    }
  } while (diff > eps_func() * eps_func());

  if (diff <= eps_func() * eps_func()){
    return true;
  } else {
    return false;
  }
}

bool SelfConsistentIntegrand2Bilayer::find_solution_nm(double& _delta, double& _mu, bool verbose){
  std::cout << "Using the Nelder-Mead method." << std::endl;    
  nm_.reset();
  nm_.set_tolerance(eps_func() * eps_func());
    
  /* Initial parameters */
  std::vector<vec> xs(n_sc_params+1);
  vec x0{_delta, _mu};
  vec x1{_delta*1.01, _mu*1.01};
  vec x2{_delta*1.01, _mu*0.99};
  xs[0] = x0;
  xs[1] = x1;
  xs[2] = x2;
  nm_.init_x(xs);
  nm_.sort();

  bool optimized = true;
  for(std::size_t t=0; t < max_iter(); ++t){
    // nm_.output(std::cout);
    nm_.step();
    if (nm_.is_terminated()) { break; }      
    if (t == max_iter() - 1){ optimized = false; }
  }
    
  double f_opt = 0.;
  nm_.get_result(x0, f_opt);
  _delta = x0[0];
  _mu = x0[1];
  if (optimized){
    return true;
  } else {
    std::cout << "Number of iteration reaches the limit " << max_iter() << std::endl;
    std::cout << "The remained error squared is " << std::setprecision(precision()) << f_opt << std::endl;
    return false;
  }
}

bool SelfConsistentIntegrand2Bilayer::find_solution_using_1d_solver(double& _delta, double& _mu, bool verbose){
  double delta0 = _delta;
  double mu0 = _mu;
  double diff1 = 0., diff2 = 0.;

  /* Precision */
  int prec = precision2();
  int pw = prec + 10;
  
  /* Testing the two alternating algorithms. */
  for(std::size_t t=0; t < max_iter(); ++t){
    /* delta -> mu */    
    _mu = calc_chemical_potential_bilayer3(L(), *ts(), filling(), T(), _delta, cbp_, continuous_k(), false);

    if (verbose) {
      std::cout << std::setw(pw) << t << std::setw(pw) << _delta << std::setw(pw) << _mu << std::endl;
    }

    /* mu -> delta */
    if (non_zero_delta()) {
      _delta = solve_self_consistent_eq_bilayer2(L(), *ts(), U(), _mu, T(), cbp_, continuous_k());
    }
    diff1 = calc_diff(_delta, _mu);
    
    if (verbose) {
      std::cout << std::setw(pw) << t << std::setw(pw) << _delta << std::setw(pw) << _mu << std::setw(pw) << diff1 << std::endl;
    }

    /* Checking if it is converged. */
    if (diff1 <= eps_func() * eps_func()) {
      return true;
    }

    /* Reached the maximum iteration step. */
    if (t == max_iter() - 1) {
      std::cout << "Number of iteration reaches the limit " << max_iter() << std::endl;
      std::cout << "The remained error squared is " << std::setprecision(precision()) << diff1 << std::endl;      
    }    
  }

  /* Storing */
  double delta1 = _delta;
  double mu1 = _mu;
  
  /* Resetting */
  _delta = delta0;
  _mu = mu0;
  
  for(std::size_t t=0; t < max_iter(); ++t){
    /* mu -> delta */
    if (non_zero_delta()) {    
      _delta = solve_self_consistent_eq_bilayer2(L(), *ts(), U(), _mu, T(), cbp_, continuous_k());
    }
    
    if (verbose){
      std::cout << std::setw(pw) << t << std::setw(pw) << _delta << std::setw(pw) << _mu << std::endl;
    }

    /* delta -> mu */    
    _mu = calc_chemical_potential_bilayer3(L(), *ts(), filling(), T(), _delta, cbp_, continuous_k(), false);    
    diff2 = calc_diff(_delta, _mu);
    
    if (verbose){
      std::cout << std::setw(pw) << t << std::setw(pw) << _delta << std::setw(pw) << _mu << std::setw(pw) << diff2 << std::endl;
    }

    /* Checking if it is converged. */
    if (diff2 <= eps_func() * eps_func()) {
      return true;
    }

    /* Reached the maximum iteration step. */
    if (t == max_iter() - 1) {
      std::cout << "Number of iteration reaches the limit " << max_iter() << std::endl;    
      std::cout << "The remained error squared is " << std::setprecision(precision()) << diff2 << std::endl;
    }    
  }  

  /* Storing */
  double delta2 = _delta;
  double mu2 = _mu;

  /* Choosing one solution. */
  if (diff1 < diff2) {
    _delta = delta1;
    _mu = mu1;
  } else {
    _delta = delta2;
    _mu = mu2;
  }
  
  return true;
}

bool SelfConsistentIntegrand2Bilayer::find_solution_default(double& _delta, double& _mu, bool verbose){
  /* Precision */
  int prec = precision2();
  int pw = prec + 10;
  
  for(std::size_t t=0; t < max_iter(); ++t){
    double delta2 = _delta;
    double mu2 = _mu;  
    /* delta -> mu */
    _mu = calc_chemical_potential_bilayer3(L(), *ts(), filling(), T(), delta2, cbp_, continuous_k(), false);

    if (non_zero_delta()) {
      /* mu -> delta */      
      _delta = solve_self_consistent_eq_bilayer2(L(), *ts(), U(), mu2, T(), cbp_, continuous_k());
    }

    /* Difference */
    double diff = calc_diff(_delta, _mu);
    
    if (verbose) {
      std::cout << std::setw(pw) << t << std::setw(pw) << _delta << std::setw(pw) << _mu << std::setw(pw) << diff << std::endl;
    }

    /* Checking if it is converged. */
    if (diff <= eps_func() * eps_func()) {
      return true;
    }

    /* Reached the maximum iteration step. */
    if (t == max_iter() - 1) {
      std::cout << "Number of iteration reaches the limit " << max_iter() << std::endl;
      std::cout << "The remained error squared is " << std::setprecision(precision()) << diff << std::endl;      
    }    
  }

  return true;
}

bool SelfConsistentIntegrand2Bilayer::find_solution(double& _delta, double& _mu, bool verbose){
  if (half_filling() && T_equal_to_0()){
    /* Checking if the charge gap of the noninteracting system is finite. */
    double delta0 = 0.;
    double ch_gap, mu0;
    std::tie(ch_gap, mu0) = calc_charge_gap_bilayer(L(), *ts(), delta0);

    /* Gapped case */
    if (ch_gap > 0.) {
      /* delta for mu0 should be correct. */
      _delta = solve_self_consistent_eq_bilayer2(L(), *ts(), U(), mu0, T(), cbp_, continuous_k());

      /* Updating mu */
      std::tie(ch_gap, _mu) = calc_charge_gap_bilayer(L(), *ts(), _delta);
      return true;
    }    
  }

  /* Two-dimensional root-finding. */
  if (use_gsl()){
    return find_solution_gsl(_delta, _mu, verbose);
  } else if (use_1d_solver()) {
    return find_solution_using_1d_solver(_delta, _mu, verbose);
  } else if (use_NewtonRaphson()) {
    return find_solution_nr(_delta, _mu, verbose);
  } else if (use_NelderMead()) {
    return find_solution_nm(_delta, _mu, verbose);
  } else {
    return find_solution_default(_delta, _mu, verbose);
  }
}

double SelfConsistentIntegrand2Bilayer::calc() const {
  double sum = 0;
  if ( continuous_k() ) {
    /* For Cuba */
    double epsrel = eps_func();
    double epsabs = eps_func();
    int ncomp = 1;    
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp_.NDIM, ncomp, sc_integrand_bilayer2, cbp_.userdata, cbp_.nvec, epsrel, epsabs, cbp_.flags, cbp_.mineval, cbp_.maxeval, cbp_.key, cbp_.statefile, cbp_.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    // for(int comp = 0; comp < ncomp; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);      

    sum = integral[0] / 4.;
  } else { /* Integral for a finite-size system */
    double k1 = 2. * M_PI / (double)L();
  
    for(int z=0; z < 2; ++z){    
      double kz = M_PI * z;
      for(int x=-L()/2; x < L()/2; ++x){
	double kx = k1 * x;
	for(int y=-L()/2; y < L()/2; ++y){
	  double ky = k1 * y;      

	  /* Checking if the wavevector is inside the BZ. */
	  double factor = BZ_factor_square_half_filling(kx, ky);
	  if ( std::abs(factor) < 1e-12 ) { continue; }

	  double qvec[3] = {kx, ky, kz};
	  sum += factor * (this->*integrand_ptr)(qvec);
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L() * L() * 2;  
    sum /= (double)(n_sites);    
  }
  
  return sum;
}

double self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double delta){
  double eps = 1e-12;
  
  /* Monotonically decreasing as a function of delta */
  double k1 = 2. * M_PI / (double)L;
  double sum = 0;

  for(int z=0; z < 2; ++z){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; ++x){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; ++y){
	double ky = k1 * y;      
	double e_eps = 1e-12;

	/* Checking if the wavevector is inside the BZ. */
	double factor = BZ_factor_square_half_filling(kx, ky);
	if ( std::abs(factor) < 1e-12 ) { continue; }

	/* Sum of the Fourier transformed hoppings */
	double ek1 = ts.ek1(kx, ky, kz);
	double ek2 = ts.ek2(kx, ky, kz);
	double ek3 = ts.ek3(kx, ky, kz);

	cx_double ak_up = calc_ak_up_in_minus(ek1, ek2, ek3, delta);
	cx_double ak_down = calc_ak_down_in_minus(ek1, ek2, ek3, delta);

	sum += factor * (std::norm(ak_up) - std::norm(ak_down));
      } /* end for y */
    } /* end for x */
  } /* end for z */

  int n_sites = L * L * 2;  
  sum /= (double)( n_sites * delta );

  /* The self-consistent condition: sum = 1/U */
  return sum;
}

double solve_self_consistent_eq_bilayer(int L, hoppings_bilayer const& ts, double U){
  std::cout << "Finding a self-consistent solution for U=" << U << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_bilayer, L, std::ref(ts), _1 );

  BinarySearch bs;
  bool sol_found = bs.find_solution( delta, target, scc );
  
  if ( sol_found ) {
    return delta;
  } else {
    return 0;
  }
}

double self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double mu, double T, double delta, CubaParam const& cbp, bool continuous_k){
  /* Monotonically decreasing as a function of delta */
  double sum = 0;

  /* Setting the parameters. */
  scib.set_parameters(ts, mu, delta, T);
    
  if ( continuous_k ) {
    /* For Cuba */
    double epsrel = 1e-10;
    double epsabs = 1e-10;    
    int ncomp = 1;
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, sc_integrand_bilayer, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);    
    // for(int comp = 0; comp < cbp.NCOMP; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);

    sum = integral[0] / ( 4. * delta );
  } else {
    double eps = 1e-12;
    double k1 = 2. * M_PI / (double)L;
    
    for(int z=0; z < 2; ++z){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; ++x){
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; ++y){
	  double ky = k1 * y;      

	  /* Checking if the wavevector is inside the BZ. */
	  double factor = BZ_factor_square_half_filling(kx, ky);
	  if ( std::abs(factor) < 1e-12 ) { continue; }

	  /* Sum of the Fourier transformed hoppings */
	  cx_double ek1 = ts.ek1(kx, ky, kz);
	  double zki_over_delta = zk_over_delta(ek1, ts.tz, kz, delta);

	  /* Fermi density */
	  double Em, Ep;
	  std::tie(Em, Ep) = calc_single_particle_energy2(ts, kx, ky, kz, delta);    
	  double mu = 0.5 * ( Em + Ep );
	  double n_minus = fermi_density(Em, kB*T, mu);
	  double n_plus = fermi_density(Ep, kB*T, mu);	    
	  sum += factor * zki_over_delta * ( n_minus - n_plus );	  
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;
    sum /= (double)( n_sites );    
  }  

  /* The self-consistent condition: sum = 1/U */
  return sum;
}

double solve_self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double U, double mu, double T, CubaParam const& cbp, bool continuous_k){
  std::cout << "Solving the self-consistent equation for U = " << U << ", mu = " << mu << ", T = " << T << std::endl;  
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind(self_consistent_eq_bilayer2, L, std::ref(ts), mu, T, _1, std::ref(cbp), continuous_k);

  BinarySearch bs(continuous_k);
  bool sol_found = bs.find_solution(delta, target, scc);

  if (sol_found) {
    return delta;
  } else {
    return 0;
  }
}

std::tuple<double, double> solve_self_consistent_eqs_bilayer(rpa::parameters const& pr, int L, hoppings_bilayer2 const& ts, double U, double filling, double T, bool continuous_k, double delta_i, double mu_i){
  /* Checking if it is in the ordered phase. */
  bool non_zero_delta;
  if (T < 1e-15) {
    std::cout << "Calculating the critical U." << std::endl;
    double Uc = find_critical_U_bilayer(pr);
    std::cout << "Uc = " << Uc << std::endl;
    
    if (U > Uc) {
      non_zero_delta = true;
    } else {
      non_zero_delta = false;
    }
  } else {
    std::cout << "Calculating the critical T." << std::endl;
    double Tc = find_critical_T_bilayer(pr);
    std::cout << "Tc = " << Tc << std::endl;
    T = pr.calc_T(Tc);
    std::cout << "T = " << T << std::endl;
    
    if (T < Tc) {
      non_zero_delta = true;
    } else {
      non_zero_delta = false;
    }
  }

  /* Disordered state */
  double delta1, mu1;
  std::tie(delta1, mu1) = solve_self_consistent_eqs_bilayer2(pr, L, ts, U, filling, T, continuous_k, false, delta_i, mu_i);
  if (non_zero_delta) {
    /* Comparing the free energies. */
    sci2b.set_input(delta1, mu1);    
    double f1 = sci2b.calc_energy();
    std::cout << "Free energy of the disordered state = " << f1 << std::endl;
    
    /* Ordered state */    
    double delta2, mu2;
    std::tie(delta2, mu2) = solve_self_consistent_eqs_bilayer2(pr, L, ts, U, filling, T, continuous_k, true, delta_i, mu_i);
    sci2b.set_input(delta2, mu2);
    double f2 = sci2b.calc_energy();
    std::cout << "Free energy of the ordered state = " << f2 << std::endl;
    
    if (f1 < f2) {
      return std::make_tuple(delta1, mu1);
    } else {
      return std::make_tuple(delta2, mu2);
    }
  } else {
    return std::make_tuple(delta1, mu1);
  }
}

std::tuple<double, double> solve_self_consistent_eqs_bilayer2(rpa::parameters const& pr, int L, hoppings_bilayer2 const& ts, double U, double filling, double T, bool continuous_k, bool non_zero_delta, double delta_i, double mu_i){
  std::cout << "Finding a self-consistent solution for U = " << U << ", T = " << T << std::endl;
  
  /* Initial values */
  double delta = 0.;
  if (non_zero_delta) {
    delta = delta_i;  
  } else {
    delta = 0.;
  }
  double mu = mu_i;
  
  /* Set parameters */
  sci2b.set_parameters(pr, L, ts, U, filling, T, delta, mu, continuous_k, non_zero_delta);

  /* Find a solution */
  bool verbose = true;
  // bool verbose = false;  
  
  bool sol_found = sci2b.find_solution(delta, mu, verbose);  

  if ( !sol_found ) {
    std::cerr << "A self-consistent solution was not found within the tolerance." << std::endl;
  }
  
  return std::make_tuple(delta, mu);  
}

void solve_self_consistent_eqs_bilayer_T(path& base_dir, rpa::parameters const& pr){  
  /* Getting parameters */
  int L = pr.L;
  double U = pr.U;
  double filling = pr.filling;
  bool continuous_k = pr.continuous_k;
  
  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);
    
  /* Ts */
  double T_min = pr.T_min;
  double T_max = pr.T_max;
  double T_delta = pr.T_delta;  
  
  int n_Ts = int((T_max - T_min)/T_delta+0.5) + 1;
  std::vector<double> Ts(n_Ts);
  for(int o=0; o < n_Ts; ++o){ Ts[o] = T_min + T_delta * o; }

  /* Output */
  ofstream out_sc;;
  out_sc.open(base_dir / "self_consistent_solution-T.text");
  out_sc << "# T   delta   mu" << std::endl;
  
  /* Precision */
  int prec = 15;
  int pw = prec + 10;
  out_sc << std::setprecision(prec);
  
  /* Tc */
  double Tc = find_critical_T_bilayer(pr);
  
  /* For each T */
  std::cout << "Solving the self-consistent equations for U = " << U << std::endl;  
  for(int Ti=0; Ti < n_Ts; ++Ti){
    double T = Ts[Ti];
    bool non_zero_delta = T < Tc;
    double delta, mu;
    std::tie(delta, mu) = solve_self_consistent_eqs_bilayer( pr, L, *ts, U, filling, T, continuous_k, non_zero_delta );

    /* Output */
    out_sc << T << std::setw(pw) << delta << std::setw(pw) << mu << std::endl;
  }
  
  out_sc.close();
}

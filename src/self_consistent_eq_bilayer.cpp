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
#include "find_critical_T.h"

/* Creating an instance */
SelfConsistentIntegrandBilayer scib;
SelfConsistentIntegrand2Bilayer sci2b;

/* Wrapper functions */
int sc_integrand_bilayer(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return scib.calc(ndim, xx, ncomp, ff, userdata);
}

int sc_integrand_bilayer2(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  sci2b.integral(ndim, xx, ncomp, ff, userdata);
  return 0;
}

/* Member functions of SelfConsistentIntegrand2Bilayer */
SelfConsistentIntegrand2Bilayer::SelfConsistentIntegrand2Bilayer():SelfConsistentIntegrand2(&hb_),nr_(n_sc_params){}

void SelfConsistentIntegrand2Bilayer::set_parameters(rpa::parameters const& pr, int _L, hoppings_bilayer2 const& _hb, double _U, double _filling, double _T, double _delta, double _mu, bool _continuous_k, bool _non_zero_delta){
  SelfConsistentIntegrand2::set_parameters(_L, _U, _filling, _T, _delta, _mu, _continuous_k, _non_zero_delta);
  SelfConsistentIntegrand2::set_eps_func(pr.epsfunc);
  hb_ = _hb;
  cbp_.set_parameters(pr);
  nr_.mod_prefactor = pr.mod_prefactor;
}

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

double SelfConsistentIntegrand2Bilayer::calc_diff(){
  if ( non_zero_delta() ) {
    return std::norm(nr_.F(0)) + std::norm(nr_.F(1));
  } else {
    return std::norm(nr_.F(1));
  }
}

void SelfConsistentIntegrand2Bilayer::update_parameters(int64_t niter, double& _delta, double& _mu){
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

bool SelfConsistentIntegrand2Bilayer::find_solution(double& delta, double& mu, bool verbose){
  double delta2 = delta;
  double mu2 = mu;
  int64_t niter = 0;
  double diff = 0;
  do {
    ++niter;
    if ( niter == max_iter() ) {
      std::cerr << "Number of iteration reaches the limit " << max_iter() << std::endl;
      std::cerr << "The remained error squared is " << std::to_string(diff) << std::endl;      
      return false;
    }
    delta = delta2;
    mu = mu2;
    set_input(delta, mu);    
    update_parameters(niter, delta2, mu2);
    diff = calc_diff();
    if ( verbose ) {
      if ( non_zero_delta() ) {
	std::cout << niter << "  "
		  << delta << "  " << delta2 << "  " << mu << "  " << mu2 << "  "
		  << diff << "  "
		  << nr_.F(0) << "  " << nr_.F(1) << "  "
		  << nr_.J(0,0) << "  " << nr_.J(0,1) << "  " << nr_.J(1,0) << "  " << nr_.J(1,1) << std::endl;
      } else {
	std::cout << niter << "  " << delta << "  "
		  << mu << "  " << mu2 << "  "
		  << diff << "  "
	          << nr_.F(1) << "  "
		  << nr_.J(1,1)
		  << std::endl;
      }
    }
  } while ( diff > eps_func() * eps_func() );
  return true;  
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
	
	  /* Checking if k is inside/outside the BZ. */
	  double mu_free = 0;  /* Assume at half filling */
	  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
	  if ( e_free > mu_free + eps() ) continue;

	  /* Prefactor */
	  double factor = 1.;
	
	  /* On the zone boundary */
	  if ( std::abs(e_free - mu_free) < eps() ) {	
	    factor = 0.5;
	  }

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
	
	/* Checking if k is inside/outside the BZ. */
	double mu_free = 0;  /* Assume at half filling */
	double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
	if ( e_free > mu_free + eps ) continue;

	/* Prefactor */
	double factor = 1.;
	
	/* On the zone boundary */
	if ( std::abs(e_free - mu_free) < eps ) {	
	  factor = 0.5;
	}

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

double self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double filling, double T, double delta, CubaParam const& cbp, bool continuous_k){  
  /* Monotonically decreasing as a function of delta */
  double sum = 0;

  /* Changing the parameters */
  double mu = calc_chemical_potential_bilayer3(L, ts, filling, T, delta, cbp, continuous_k, false);
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
	
	  /* Checking if k is inside/outside the BZ. */
	  double mu_free = 0;  /* Assume at half filling */
	  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
	  if ( e_free > mu_free + eps ) continue;

	  /* Prefactor */
	  double factor = 1.;
	
	  /* On the zone boundary */
	  if ( std::abs(e_free - mu_free) < eps ) {	
	    factor = 0.5;
	  }

	  /* Sum of the Fourier transformed hoppings */
	  cx_double ek1 = ts.ek1(kx, ky, kz);
	  double zki_over_delta = zk_over_delta(ek1, ts.tz, kz, delta);

	  /* Fermi density */
	  double Em, Ep;
	  std::tie(Em, Ep) = calc_single_particle_energy2(ts, kx, ky, kz, delta);    
	  double mu = 0.5 * ( Em + Ep );
	  double n_minus = fermi_density(Em, kB*T, mu);
	  double n_plus = fermi_density(Ep, kB*T, mu);	    

	  // // for check
	  // double zki = zk(ek1, ts.tz, kz, delta);	  
	  // sum += factor * zki;
	  
	  sum += factor * zki_over_delta * ( n_minus - n_plus );	  
	  
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    // sum /= (double)( n_sites * delta );
    sum /= (double)( n_sites );    
  }  

  /* The self-consistent condition: sum = 1/U */
  return sum;
}

double solve_self_consistent_eq_bilayer2(int L, hoppings_bilayer2 const& ts, double U, double filling, double T, CubaParam const& cbp, bool continuous_k){
  
  std::cout << "Finding a self-consistent solution for U = " << U << ", T = " << T << std::endl;
  double target = 1. / U;
  double delta = 0.45 * U;

  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_bilayer2, L, std::ref(ts), filling, T, _1, std::ref(cbp), continuous_k );

  BinarySearch bs(continuous_k);
  bool sol_found = bs.find_solution( delta, target, scc );

  if ( sol_found ) {
    return delta;
  } else {
    return 0;
  }  
}

std::tuple<double, double> solve_self_consistent_eqs_bilayer(rpa::parameters const& pr, int L, hoppings_bilayer2 const& ts, double U, double filling, double T, bool continuous_k, bool non_zero_delta, double delta_i = 0.01, double mu_i = 0.){  
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

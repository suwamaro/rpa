/*****************************************************************************
*
* Functions for calculating the wave function
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <complex>
#ifdef WITH_OpenMP
#include <omp.h>
#endif
#include "calc_intensity.h"
#include "calc_spectrum.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"
#include "calc_chemical_potential.h"
#include "calc_wave_func.h"
#include "BinarySearch.h"

/* Creating an instance */
PhiDerIntegrandBilayer pdib;
WaveFuncIntegrandBilayer wfib;

/* Member functions of PhiDerIntegrand */
double PhiDerIntegrand::integrand(double omega, double delta, double zk, cx_double bk){
  if ( delta == 0 ) {
    /* For U <= Uc */
    double b = std::abs(bk);
    double g2 = std::pow(omega / (2.*b), 2);
    return 1. / ( std::pow(b,3) * std::pow( 1. - g2, 2) );
  } else {
    /* For U > Uc */

    // // for check
    // std::cout << "Psider: " << zk << "  " << delta << "  " << omega << std::endl;
    
    return std::pow(zk,3) * (1. - zk*zk) / std::pow(1 - std::pow(omega*zk/(2.*delta), 2), 2);
  }
}

/* Member functions of PhiDerIntegrandBilayer */
void PhiDerIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double omega, double delta){
  hb_ = h;
  omega_ = omega;
  delta_ = delta;
}


int PhiDerIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
  /* Reset */
  ff[0] = 0;
  
  /* Wavenumbers */
  double k1 = xx[0] * 2 * M_PI;
  double k2 = xx[1] * 2 * M_PI;
  
  double kx = 0.5 * (k2 + k1);
  double ky = 0.5 * (k2 - k1);
  
  /* Sum over kz */
  for(int z=0; z < 2; z++){       
    double kz = M_PI * z;	  
    
    cx_double ek1 = ts()->ek1(kx, ky, kz);
    double zki = zk(ek1, ts()->tz, kz, delta());
    cx_double bki = bk(up_spin, ek1, ts()->tz, kz);  // spin does not matter.
    ff[0] += integrand(omega(), delta(), zki, bki);
  }
  
  return 0;   
}

/* Member functions of WaveFuncIntegrandBilayer */
cx_double WaveFuncIntegrandBilayer::integrand(cx_double xk, double zk, cx_double bk, double ek_plus, double ek_minus, cx_double phase) const {
  return integrand(xk, zk, bk, ek_plus, ek_minus, phase, sublattice());
}

cx_double WaveFuncIntegrandBilayer::integrand(cx_double xk, double zk, cx_double bk, double ek_plus, double ek_minus, cx_double phase, int sublattice) const {
  if ( sublattice == 0 ) {
    if ( largeUlimit() ) {
      return 0;
    } else {
      double ek_diff = ek_plus - ek_minus;
      return (double)spin() * (1. - zk*zk) / (2.*sqrt(psider())) * ek_diff / (omega()*omega() - ek_diff*ek_diff) * phase;
    }
  } else if ( sublattice == 1 ) {
    if ( largeUlimit() ) {
      return 2. * std::conj(bk) / sqrt(psider()) / ( scaling_prefactor() + 2. * (std::norm(bk) - min_bk_sq()) ) * phase;        
    } else {
      double ek_diff = ek_plus - ek_minus;
      return (double)spin() * std::conj(xk) * sqrt(1. - zk*zk) / (2.* sqrt(psider())) * ( (1.-spin()*zk) / (omega() - ek_diff) + (1.+spin()*zk) / (omega() + ek_diff) ) * phase;
    }
  } else {
    std::cerr << "\"sublattice_\" has to be 0 or 1." << std::endl;
    std::exit(EXIT_FAILURE);
  }  
}

void WaveFuncIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, bool largeUlimit, double scaling_prefactor, int spin, double omega, double psider, double delta, int *diff_r, int sublattice, double min_bk_sq){
  hb_ = h;
  largeUlimit_ = largeUlimit;
  scaling_prefactor_ = scaling_prefactor;
  spin_ = spin;
  omega_ = omega;
  psider_ = psider;
  delta_ = delta;
  for(int d=0; d < 3; d++){ diff_r_[d] = diff_r[d]; }
  sublattice_ = sublattice;
  min_bk_sq_ = min_bk_sq;
}


int WaveFuncIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
  cx_double sum = 0;
  
  /* Wavenumbers */
  double k1 = xx[0] * 2 * M_PI;
  double k2 = xx[1] * 2 * M_PI;
  
  double kx = 0.5 * (k2 + k1);
  double ky = 0.5 * (k2 - k1);
  
  /* Sum over kz */
  for(int z=0; z < 2; z++){       
    double kz = M_PI * z;	  
    
    cx_double ek1 = ts()->ek1(kx, ky, kz);
    cx_double ek23 = ts()->ek23(kx, ky, kz);
    cx_double ekz = ts()->ekz(kx, ky, kz);    
    cx_double xki = xk(spin(), ek1, ts()->tz, kz, delta());
    double zki = zk(ek1, ts()->tz, kz, delta());
    cx_double bki = bk(spin(), ek1, ts()->tz, kz);
    double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts()->tz, kz, delta());    
    double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts()->tz, kz, delta()); 
    cx_double pure_i(0,1);
    double inner_prod = kx * diff_r(0) + ky * diff_r(1) + kz * diff_r(2);
    cx_double phase = exp(pure_i*inner_prod);
    sum += integrand(xki, zki, bki, ek_plus, ek_minus, phase);
  }

  ff[0] = std::real(sum);
  ff[1] = std::imag(sum);
  
  return 0;   
}

int pd_integrand_bilayer(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return pdib.calc(ndim, xx, ncomp, ff, userdata);
}

int wf_integrand_bilayer(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return wfib.calc(ndim, xx, ncomp, ff, userdata);
}

double pole_eq_bilayer(int L, hoppings_bilayer2 const& ts, double omega, double mu, double U, double T, double delta, CubaParam const& cbp, MatElemF const& me_F, bool continuous_k, std::string const& mode){
  arma::cx_mat chi0_pm(NSUBL,NSUBL,arma::fill::zeros);
  arma::cx_mat chi0_zz(NSUBL,NSUBL,arma::fill::zeros);

  /* Calculating the bare response functions */
  std::tie(chi0_pm, chi0_zz) = calc_bare_response_bilayer(L, ts, mu, U, T, delta, cbp, me_F, omega, continuous_k);
    
  cx_double det = 0;
  if ( mode == "transverse" ) {
    det = arma::det(arma::eye<arma::cx_mat>(NSUBL,NSUBL) - U * chi0_pm);    
  } else if ( mode == "longitudinal" ) {
    det = arma::det(arma::eye<arma::cx_mat>(NSUBL,NSUBL) - 0.5 * U * chi0_zz);
  } else {
    std::cerr << "\"mode\" has to be \"transverse\" or \"longitudinal\"." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  return std::real(det);  /* Assume that the imaginary part is 0. */
}

double solve_pole_eq_bilayer(int L, hoppings_bilayer2 const& ts, double mu, double U, double T, double delta, CubaParam const& cbp, MatElemF const& me_F, bool continuous_k, std::string const& mode, double upper, bool return_upper, bool verbose){
  if ( verbose ) {
    std::cout << "Solving the pole equation for U=" << U << std::endl;
    std::cout << "Upper = " << upper << std::endl;
  }
  
  double target = 0;
  double omega = 0.5 * upper;

  using std::placeholders::_1;
  auto pe = std::bind( pole_eq_bilayer, L, std::ref(ts), _1, mu, U, T, delta, std::ref(cbp), std::ref(me_F), continuous_k, mode );
  BinarySearch bs(continuous_k);
  bs.set_x_MIN(0);
  double omega_eps = 1e-5;
  // double omega_eps = 1e-7;  
  bs.set_x_MAX(upper - omega_eps);
  bool sol_found = bs.find_solution( omega, target, pe, false, 0, verbose );
  
  if ( sol_found ) {
    return omega;
  } else {
    if ( return_upper ) {
      return upper;
    } else {
      return 0;
    }
  }
}

std::tuple<double, double, double> calc_gap_bilayer(int L, hoppings_bilayer2 const& ts, double mu, double U, double T, double delta, CubaParam const& cbp, MatElemF const& me_F, bool continuous_k, bool return_upper, bool verbose){
  std::cout << "Calculating the excitation gaps." << std::endl;
  double upper = calc_band_gap_bilayer(L, ts, delta, me_F.qx(), me_F.qy(), me_F.qz());
  double omega_T = solve_pole_eq_bilayer(L, ts, mu, U, T, delta, cbp, me_F, continuous_k, "transverse", upper, return_upper, verbose);
  double omega_L = solve_pole_eq_bilayer(L, ts, mu, U, T, delta, cbp, me_F, continuous_k, "longitudinal", upper, return_upper, verbose);
  return std::make_tuple(omega_T, omega_L, upper);
}

double calc_Psider(int L, hoppings_bilayer2 const& ts, double mu, double omega, double delta, CubaParam const& cbp, bool continuous_k){
  double Psider = 0;
  if ( continuous_k ) {
    /* Changing the parameters */
    pdib.set_parameters(ts, omega, delta);

    /* For Cuba */
    double epsrel = 1e-10;
    double epsabs = 1e-10;    
    int ncomp = 1;
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, pd_integrand_bilayer, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);    
    // for(int comp = 0; comp < ncomp; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);

    if ( delta == 0 ) {
      Psider = integral[0] * omega / std::pow(2, 3);
    } else {
      Psider = integral[0] * omega / std::pow(2*delta, 3);
    }
  } else {
    double eps = 1e-12;
    double k1 = 2. * M_PI / (double)L;
  
    for(int z=0; z < 2; z++){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; x++){
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; y++){
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
	  double zki = zk(ek1, ts.tz, kz, delta);
	  cx_double bki = bk(up_spin, ek1, ts.tz, kz);  // spin does not matter.
	  Psider += factor * PhiDerIntegrand::integrand(omega, delta, zki, bki);

	  // // for check
	  // std::cout << x << "  " << y << "  " << z << "  " << omega << "  " << delta << "  " << zki << "  " << bki << "  " << factor << "  " << Psider << std::endl;
	  
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;
    int n_unit_cell = n_sites / 2;

    if ( delta == 0 ) {
      Psider *=  2. * omega / std::pow(2., 3) / (double)n_unit_cell;
    } else {
      Psider *=  2. * omega / std::pow(2*delta, 3) / (double)n_unit_cell;
    }
  }
  
  return Psider;
}

double calc_min_bk_sq_bilayer(int L, hoppings2 const& ts){
  /* Finite size */
  double k1 = 2. * M_PI / (double)L;
  double min_bk_sq = std::numeric_limits<double>::max();

  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;
	cx_double ek1 = ts.ek1(kx, ky, kz);
	cx_double bki = bk(up_spin, ek1, ts.tz, kz);  // spin does not matter.
	min_bk_sq = std::min( min_bk_sq, std::norm(bki) );
      } /* end for y */
    } /* end for x */
  } /* end for z */

  return min_bk_sq;
}
 
double calc_weight_distribution_bilayer(int L, hoppings_bilayer2 const& ts, CubaParam const& cbp, MatElemF const& me_F, bool continuous_k, int *diff_r){    
  cx_double wavefunc = 0;  
  if ( continuous_k ) {
    /* For Cuba */
    double epsrel = 1e-10;
    double epsabs = 1e-10;    
    int ncomp = 2;
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, wf_integrand_bilayer, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);    
    // for(int comp = 0; comp < ncomp; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);

    wavefunc = { integral[0] / 2., integral[1] / 2. };
  } else {
    double eps = 1e-12;
    double k1 = 2. * M_PI / (double)L;
  
    for(int z=0; z < 2; z++){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; x++){
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; y++){
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
	  
	  cx_double ek1 = ts.ek1(kx, ky, kz);
	  cx_double ek23 = ts.ek23(kx, ky, kz);
	  cx_double ekz = ts.ekz(kx, ky, kz);    
	  cx_double xki = xk(wfib.spin(), ek1, ts.tz, kz, wfib.delta());
	  double zki = zk(ek1, ts.tz, kz, wfib.delta());
	  cx_double bki = bk(wfib.spin(), ek1, ts.tz, kz);
	  double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts.tz, kz, wfib.delta());
	  double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts.tz, kz, wfib.delta());
	  cx_double pure_i(0,1);
	  double inner_prod = kx * diff_r[0] + ky * diff_r[1] + kz * diff_r[2];
	  cx_double phase = exp(pure_i*inner_prod);
	  wavefunc += factor * wfib.integrand(xki, zki, bki, ek_plus, ek_minus, phase);
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;
    int n_unit_cell = n_sites / 2;
    wavefunc /=  (double)n_unit_cell;
  }

  cx_double pure_i(0,1);
  double inner_prod2 = me_F.qx() * diff_r[0] + me_F.qy() * diff_r[1] + me_F.qz() * diff_r[2];  // Assume that the origin is (0,0,0).
  cx_double phase2 = exp(0.5*pure_i*inner_prod2);
  wavefunc *= phase2;
  
  return std::norm(wavefunc);
}

void calc_wave_func_bilayer(path& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  int Lk = pr.Lk;
  double U = pr.U;
  double filling = pr.filling;
  double T = pr.T;
  bool continuous_k = pr.continuous_k;

  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);
  
  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  /* Calculate the chemical potential and the charge gap. */
  double delta = solve_self_consistent_eq_bilayer2( L, *ts, U, filling, T, cbp, continuous_k );  
  std::cout << "delta = " << delta << std::endl;  
  /* Assume that mu does not depend on L for integral over continuous k. */
  double ch_gap, mu;
  std::tie(ch_gap, mu) = calc_charge_gap_bilayer( L, *ts, delta );  /* Finite size */  
  // double mu = calc_chemical_potential_bilayer2( base_dir, L, *ts, delta );  /* Finite size */

  /* Output */
  ofstream out_gap;
  out_gap.open(base_dir / "gap.text");
  out_gap << "Charge gap = " << ch_gap << std::endl;
  out_gap << "mu = " << mu << std::endl;
  
  /* Precision */
  int prec = 15;

  /* MatElemF */
  MatElemF me_F(L, L, 2, 1, NSUBL);
  
  /* Setting the ordering vector */
  me_F.set_q( M_PI, M_PI, M_PI );
  if ( !continuous_k ) {
    /* The polarizations are calculated in advance. */
    me_F.set_table( *ts, delta );
  }

  /* Calculating gaps */
  double omega_T = 0, omega_L = 0, omega_ph = 0;
  std::tie(omega_T, omega_L, omega_ph) = calc_gap_bilayer(L, *ts, mu, U, T, delta, cbp, me_F, continuous_k);

  /* Output */
  out_gap << "omega_T = " << omega_T << std::endl;
  out_gap << "omega_L = " << omega_L << std::endl;
  out_gap << "omega_ph = " << omega_ph << std::endl;
  
  /* Calculating Psi'*/
  double Psider = calc_Psider(L, *ts, mu, omega_L, delta, cbp, continuous_k);

  /* Output */
  out_gap << "Psider = " << Psider << std::endl;

  /* Output */
  ofstream out_prob;
  out_prob.open(base_dir / "weight-distribution.text");
  out_prob << "# x    y    z    P(r)" << std::endl;

  /* Calculating the minimum bk square */
  double min_bk_sq = calc_min_bk_sq_bilayer(L, *ts);
  
  /* Calculating the weight distribution of the exciton eigenstate */
  int spin = -1; // down electron and down hole
  int dL = 10;
  for(int dx=-dL; dx <= dL; dx++){
    for(int dy=-dL; dy <= dL; dy++){
      for(int dz=0; dz <= 1; dz++){
	int diff_r[] = {dx, dy, dz};
	int sublattice = (dx+dy+dz) & 1;
	
	/* Setting the parameters */	  
	wfib.set_parameters(*ts, pr.largeUlimit, pr.largeU_scaling_prefactor, spin, omega_L, Psider, delta, diff_r, sublattice, min_bk_sq);

	/* Calculating the weight distribution */
	double prob = calc_weight_distribution_bilayer(L, *ts, cbp, me_F, continuous_k, diff_r);

	/* Output */
	out_prob << dx << std::setw(10) << dy << std::setw(10) << dz << std::setw(20) << prob << std::endl;
      }
    }
  }

  out_gap.close();
  out_prob.close();
}

void check_wave_func_bilayer(path& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  int Lk = pr.Lk;
  double U = pr.U;
  double filling = pr.filling;
  double T = pr.T;
  bool continuous_k = pr.continuous_k;

  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);
  
  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  /* Calculate the chemical potential and the charge gap. */
  double delta0 = 0.;
  double mu0 = calc_chemical_potential_bilayer3(L, *ts, filling, T, delta0, cbp, continuous_k, false);  
  double delta = solve_self_consistent_eq_bilayer2( L, *ts, U, mu0, T, cbp, continuous_k );  
  std::cout << "delta = " << delta << std::endl;  
  /* Assume that mu does not depend on L for integral over continuous k. */
  double ch_gap, mu;
  std::tie(ch_gap, mu) = calc_charge_gap_bilayer( L, *ts, delta );  /* Finite size */  
  // double mu = calc_chemical_potential_bilayer2( base_dir, L, *ts, delta );  /* Finite size */

  /* Output */
  ofstream out_gap;
  out_gap.open(base_dir / "gap.text");
  out_gap << "Charge gap = " << ch_gap << std::endl;
  out_gap << "mu = " << mu << std::endl;
  
  /* MatElemF */
  MatElemF me_F(L, L, 2, 1, NSUBL);
  
  /* Setting the ordering vector */
  me_F.set_q( M_PI, M_PI, M_PI );
  if ( !continuous_k ) {
    /* The polarizations are calculated in advance. */
    me_F.set_table( *ts, delta );
  }

  /* Calculating gaps */
  double omega_T = 0, omega_L = 0, omega_ph = 0;
  std::tie(omega_T, omega_L, omega_ph) = calc_gap_bilayer(L, *ts, mu, U, T, delta, cbp, me_F, continuous_k);

  /* Output */
  out_gap << "omega_T = " << omega_T << std::endl;
  out_gap << "omega_L = " << omega_L << std::endl;
  out_gap << "omega_ph = " << omega_ph << std::endl;
  
  /* Calculating Psi'*/
  double Psider = 1.0;
  // double Psider = calc_Psider(L, *ts, mu, omega_L, delta, cbp, continuous_k);

  /* Output */
  out_gap << "Psider = " << Psider << std::endl;

  /* Output */
  ofstream out_wf;
  out_wf.open(base_dir / "exciton_wave_function-k.text");
  out_wf << "# x    y    z    P(r)" << std::endl;

  /* Calculating the minimum bk square */
  double min_bk_sq = calc_min_bk_sq_bilayer(L, *ts);
  
  /* Calculating the wave function */
  int spin = -1; // down electron and down hole
  int sublattice = 0;
  int diff_r[] = {0,0,0};
  
  /* Setting the parameters */	  
  wfib.set_parameters(*ts, pr.largeUlimit, pr.largeU_scaling_prefactor, spin, omega_L, Psider, delta, diff_r, sublattice, min_bk_sq);

  double k1 = 2. * M_PI / (double)L;
  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;      
	cx_double ek1 = ts->ek1(kx, ky, kz);
	cx_double ek23 = ts->ek23(kx, ky, kz);
	cx_double ekz = ts->ekz(kx, ky, kz);    
	cx_double xki = xk(wfib.spin(), ek1, ts->tz, kz, wfib.delta());
	double zki = zk(ek1, ts->tz, kz, wfib.delta());
	cx_double bki = bk(wfib.spin(), ek1, ts->tz, kz);
	double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts->tz, kz, wfib.delta());
	double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts->tz, kz, wfib.delta());
	double wavefunc = std::real(wfib.integrand(xki, zki, bki, ek_plus, ek_minus, 1.));
	
	out_wf << kx << std::setw(15) << ky << std::setw(15) << kz << std::setw(15) << wavefunc << std::endl;
      } /* end for y */
    } /* end for x */
  } /* end for z */

  out_gap.close();
  out_wf.close();
}
			 

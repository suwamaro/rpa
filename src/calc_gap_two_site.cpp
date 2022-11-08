/*****************************************************************************
*
* Functions for calculating the gap.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_gap.h"
#include "rpa_util.h"
#include "BinarySearch.h"

cx_double calc_prefactor_bare_res_func_two_site(int sg1, int sg2, hoppings2 const& ts, double T, double k, double q, cx_double omega, double delta, double mu){
  /* k */
  cx_double ek1_k = ts.ek1(k, 0., 0.);
  double Ek = eigenenergy_HF(sg1, ek1_k, 0, 0, ts.tz, 0, delta);
  double nk = fermi_density(Ek, kB*T, mu);
  
  /* k + q */
  double k_q = k + q;
  cx_double ek1_k_q = ts.ek1(k_q, 0., 0.);  
  double Ek_q = eigenenergy_HF(sg2, ek1_k_q, 0, 0, ts.tz, 0., delta);
  double nk_q = fermi_density(Ek_q, kB*T, mu);

  /* Denominator */
  cx_double diff_E = omega - (Ek_q - Ek);
  
  /* Electron density at zero temperature */  
  // int n1 = ( ( sg1 + 1 ) >> 1 ) ^ 1;  /* -1 -> 1, 1 -> 0 */
  // int n2 = n1 ^ 1;
  // double n_diff = (double)(n2 - n1);
  double n_diff = nk_q - nk;
  
  cx_double prefactor = n_diff / diff_E;

  // for check
  std::cout << "Prefactor: " << Ek << "  " << nk << "  " << Ek_q << "  " << nk_q << "  " << diff_E << "  " << n_diff << "  " << prefactor << std::endl;
  
  return prefactor;
}

void add_to_sus_mat_two_site(hoppings2 const& ts, double T, double mu, arma::cx_mat& chi_pm, arma::cx_mat& chi_zz_u, double k, Polarization const& Pz, double delta, cx_double omega){  
  for(int sg1=-1; sg1<=1; sg1+=2){
    for(int sg2=-1; sg2<=1; sg2+=2){
      cx_double prefactor = calc_prefactor_bare_res_func_two_site(sg1, sg2, ts, T, k, Pz.qz(), omega, delta, mu);
      int sg1i = (sg1+1) >> 1;    
      int sg2i = (sg2+1) >> 1;    
      cx_double Ppm[NSUBL*NSUBL];
      cx_double Pzz[NSUBL*NSUBL];    

      if ( Pz.is_table_set() ) {
	Pz.get_Ppm(0, 0, k, sg1i, sg2i, Ppm);
	Pz.get_Pzz(0, 0, k, sg1i, sg2i, Pzz);
      } else {
	Pz.calc_polarization(ts, delta, 0, 0, k, sg1, sg2, Ppm, Pzz);
      }
      
      chi_pm(0,0) += prefactor * Ppm[0];
      chi_pm(0,1) += prefactor * Ppm[1];
      chi_pm(1,0) += prefactor * Ppm[2];
      chi_pm(1,1) += prefactor * Ppm[3];

      // for check
      std::cout << sg1 << "  " << sg2 << "  " << prefactor << "  " << Pzz[0] << "  " << Pzz[1] << "  " << Pzz[2] << "  " << Pzz[3] << std::endl;
      
      chi_zz_u(0,0) += prefactor * Pzz[0];
      chi_zz_u(0,1) += prefactor * Pzz[1];
      chi_zz_u(1,0) += prefactor * Pzz[2];
      chi_zz_u(1,1) += prefactor * Pzz[3];
    }
  }
}

std::tuple<arma::cx_mat, arma::cx_mat> calc_bare_response_two_site(hoppings_two_site const& ts, double mu, double U, double T, double delta, Polarization const& Pz, cx_double omega){
  arma::cx_mat chi0_pm(NSUBL,NSUBL,arma::fill::zeros);
  arma::cx_mat chi0_zz_u(NSUBL,NSUBL,arma::fill::zeros);
  
  add_to_sus_mat_two_site( ts, T, mu, chi0_pm, chi0_zz_u, 0, Pz, delta, omega );	  
  
  int n_units = 1;  // Number of unit cells
  chi0_pm /= (double)(n_units);
  chi0_zz_u /= (double)(n_units);
  
  /* Adding the contribution from the down spin */
  arma::cx_mat chi0_zz_d(chi0_zz_u);
  /* Assume NSUBL == 2. */  
  if ( NSUBL == 2 ) {
    chi0_zz_d.swap_rows( 0, 1 );
    chi0_zz_d.swap_cols( 0, 1 );
  } else {
    std::cerr << "NSUBL is not 2.\n";
    std::exit(EXIT_FAILURE);
  }
  arma::cx_mat chi0_zz = chi0_zz_u + chi0_zz_d;

  // for check
  std::cout << "omega = " << omega << std::endl;
  std::cout << chi0_zz_u << std::endl;
  std::cout << chi0_zz << std::endl;  
  
  return std::make_tuple(chi0_pm, chi0_zz);
}

double pole_eq_two_site(hoppings_two_site const& ts, double omega, double mu, double U, double T, double delta, Polarization const& Pz, std::string const& mode){
  arma::cx_mat chi0_pm(NSUBL,NSUBL,arma::fill::zeros);
  arma::cx_mat chi0_zz(NSUBL,NSUBL,arma::fill::zeros);

  /* Calculating the bare response functions */
  std::tie(chi0_pm, chi0_zz) = calc_bare_response_two_site(ts, mu, U, T, delta, Pz, omega);

  // for check
  std::cout << "omega = " << omega << std::endl;
  std::cout << chi0_zz << std::endl;
  
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

double solve_pole_eq_two_site(hoppings_two_site const& ts, double mu, double U, double T, double delta, Polarization const& Pz, std::string const& mode, double upper, bool return_upper, bool verbose){
  if ( verbose ) {
    std::cout << "Solving the pole equation for U=" << U << std::endl;
    std::cout << "Upper = " << upper << std::endl;
  }
  
  double target = 0;
  double omega = 0.5 * upper;

  using std::placeholders::_1;
  auto pe = std::bind( pole_eq_two_site, std::ref(ts), _1, mu, U, T, delta, std::ref(Pz), mode );
  bool continuous_k = false;
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

double calc_band_gap_two_site(hoppings_two_site const& ts, double delta, double k){
  double E_diff_min = std::numeric_limits<double>::max();

  for(int x=0; x < 2; ++x){
    double k = M_PI * x;
    cx_double ek1 = ts.ek1(k, 0., 0.);
    double E_plus = eigenenergy_HF(1, ek1, 0, 0, ts.tz, 0., delta);
    double E_minus = eigenenergy_HF(-1, ek1, 0, 0, ts.tz, 0., delta);

    double E_diff = std::max( 0., E_plus - E_minus );
    E_diff_min = std::min( E_diff_min, E_diff );
  } /* end for z */
  
  return E_diff_min;
}

std::tuple<double, double, double> calc_gap_two_site(hoppings_two_site const& ts, double mu, double U, double T, double delta, Polarization const& Pz, bool return_upper, bool verbose){
  double upper = calc_band_gap_two_site(ts, delta, Pz.qz());
  double omega_T = solve_pole_eq_two_site(ts, mu, U, T, delta, Pz, "transverse", upper, return_upper, verbose);
  double omega_L = solve_pole_eq_two_site(ts, mu, U, T, delta, Pz, "longitudinal", upper, return_upper, verbose);
  return std::make_tuple(omega_T, omega_L, upper);
}

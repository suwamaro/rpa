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
#include "array3.h"
#include "calc_intensity.h"
#include "calc_spectrum.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"
#include "calc_chemical_potential.h"
#include "calc_wave_func.h"
#include "BinarySearch.h"
#include "calc_Raman.h"
// #include "particle_hole1.h"

/* Creating an instance */
extern WaveFuncIntegrandBilayer wfib;

typedef rpa::array3 vec3;

int inner_prod(BondDelta const& d1, BondDelta const& d2){
  return d1.x * d2.x + d1.y * d2.y + d1.z * d2.z;
}

cx_double bond_to_hopping_bilayer(hoppings2 const& ts, BondDelta const& b, int g1, int g2, int sigma){
  /* Hopping amplitude from g2 to g1 */
  if (g1 == g2) {
    if (b.x + b.y + b.z & 1) {
      return cx_double(0, 0);
    } else {
      if ( b.z == 0 ) {
	return ts.tp;
      } else {
	return ts.tzp;
      }
    }
  } else {
    if (b.x + b.y + b.z & 1) {
      if ( std::abs(b.x + b.y) == 1 ) {
	return cx_double(ts.t);
      } else if ( std::abs(b.z) == 1 ) {
	if (g2 == 0) {
	  if ( sigma == up_spin ) {
	    return ts.tz;
	  } else {
	    return std::conj(ts.tz);
	  }
	} else {
	  if ( sigma == up_spin ) {
	    return std::conj(ts.tz);
	  } else {
	    return ts.tz;
	  }
	}
      } else {
	std::cerr << "Not supported hopping case." << std::endl;
	std::exit(EXIT_FAILURE);
      }
    } else {
      return cx_double(0, 0);
    }
  }
}

cx_double velocity_U1(hoppings2 const& ts, cx_mat const& Udg, cx_mat const& U, cx_mat const& U_bar, double kx, double ky, double kz, vec3 const& photon_q, BondDelta const& e_mu, std::vector<BondDelta> const& bonds, int sign_m, int sign_n, int sigma_m, int sigma_n){
  cx_double v_tot = 0;
  for(BondDelta bond: bonds){  
    int inner_prodbond = inner_prod(e_mu, bond);
    if ( inner_prodbond == 0 ) { continue; }
    
    cx_double v = 0;  
		
    /* Checking if q == 0. */
    bool q_pi = false;
    if ( photon_q.norm2() > 1e-20 ) { q_pi = true; }   // Assume q = 0 or Q
  
    /* Phases */
    double inner_prod_k = kx * bond.x + ky * bond.y + kz * bond.z;
    cx_double phase1 = exp(1i*inner_prod_k);
    cx_double phase_delta = 0;
    if ( q_pi ) {
      double inner_prod_pi_delta = photon_q[0] * bond.x + photon_q[1] * bond.y + photon_q[2] * bond.z;
      phase_delta = exp(1i*0.5*inner_prod_pi_delta);
    } else {
      phase_delta = 1.0;
    }
  
    /* Assume there is no hopping with spin flip [U(1) symmetry]. */  
    int m_idx = sign_spin_index(sign_m, sigma_m);
    int n_idx = sign_spin_index(sign_n, sigma_n);  
  
    if ( (bond.x + bond.y) & 1 ) {
      /* Different sublattices in the plane */
      for(int _spin: {up_spin, down_spin}){
	int A_idx = sublattice_spin_index(0, _spin);
	int B_idx = sublattice_spin_index(1, _spin);
	cx_double UdmA = Udg(m_idx, A_idx);      
	cx_double UdmB = Udg(m_idx, B_idx);
	cx_double UAn = U(A_idx, n_idx);    
	cx_double UBn = U(B_idx, n_idx);
	cx_double t = bond_to_hopping_bilayer(ts, bond, 0, 1, _spin);
	if ( q_pi ) {
	  cx_double epsilon = 2.0 * std::real(- t * phase1);
	  v += - 1i * phase_delta * epsilon * (UdmA * UBn - UdmB * UAn);
	} else {
	  cx_double vel = - 2.0 * std::imag(- t * phase1);  // Negative sign for the hopping amplitude.
	  v += vel * (UdmA * UBn + UdmB * UAn);	
	}
      }
    } else if (bond.z == 1) {
      /* Vertical bonds */
      for(int _spin: {up_spin, down_spin}){
	int A_idx = sublattice_spin_index(0, _spin);
	int B_idx = sublattice_spin_index(1, _spin);
	cx_double UdmA = Udg(m_idx, A_idx);      
	cx_double UdmB = Udg(m_idx, B_idx);
	cx_double UAn = U_bar(A_idx, n_idx);
	cx_double UBn = U_bar(B_idx, n_idx);
	cx_double tz = bond_to_hopping_bilayer(ts, bond, 1, 0, _spin);  // From B to A
	cx_double phase_z = exp(1i*kz);
	cx_double v0 = - 1i * phase_z * (std::conj(tz) * UdmA * UBn + tz * UdmB * UAn);
	if ( q_pi ) {
	  v += - 1i * v0;
	} else {
	  v += v0;
	}
      }    
    } else {
      /* Same sublattices */
      for(int _spin: {up_spin, down_spin}){
	int A_idx = sublattice_spin_index(0, _spin);
	int B_idx = sublattice_spin_index(1, _spin);
	cx_double UdmA = Udg(m_idx, A_idx);      
	cx_double UdmB = Udg(m_idx, B_idx);
	cx_double UAn = U(A_idx, n_idx);    
	cx_double UBn = U(B_idx, n_idx);
	cx_double t = bond_to_hopping_bilayer(ts, bond, 0, 0, _spin);
	cx_double vel = - 2.0 * std::imag(- t * phase1);  // Negative sign for the hopping amplitude.
	if ( q_pi ) {
	  v += - phase_delta * vel * (UdmA * UAn - UdmB * UBn);
	} else {
	  v += vel * (UdmA * UAn + UdmB * UBn);
	}      
      }
    }
    v_tot += (double)inner_prodbond * v;
  }
  return v_tot;
}

cx_double calc_coef_eff_Raman(int L, hoppings2 const& ts, double delta, double kx, double ky, double kz, int sigma, std::vector<BondDelta> const& bonds, BondDelta mu, BondDelta nu, vec3 const& ki, vec3 const& kf){
  /* Eigenenergy */
  cx_double ek1 = ts.ek1(kx, ky, kz);
  cx_mat Uk = gs_HF(ek1, ts.tz, kz, delta);
  cx_mat Uk_bar = gs_HF(ek1, ts.tz, kz + M_PI, delta);   // kz + M_PI
  cx_mat Uk_dg = Uk.t();

  /* Velocity */
  cx_double coef = 0.;
  if ( mu.z == 1 && nu.z == 1 ) {
    /* zz */
    double kz_bar = kz + M_PI;
    cx_double v1 = velocity_U1(ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, kf, nu, bonds, 1, 1, sigma, sigma);
    cx_double v2 = velocity_U1(ts, Uk_dg, Uk, Uk_bar, kx, ky, kz_bar, - ki, mu, bonds, 1, -1, sigma, sigma);
    
    cx_double v3 = velocity_U1(ts, Uk_dg, Uk, Uk_bar, kx, ky, kz_bar, kf, nu, bonds, -1, -1, sigma, sigma);
    cx_double v4 = velocity_U1(ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, - ki, mu, bonds, 1, -1, sigma, sigma);
    
    coef = v1 * v2 + v3 * v4;
  } else {
    cx_double v1 = velocity_U1(ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, kf, nu, bonds, 1, 1, sigma, sigma);
    cx_double v2 = velocity_U1(ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, kf, nu, bonds, -1, -1, sigma, sigma);  
    cx_double v3 = velocity_U1(ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, - ki, mu, bonds, 1, -1, sigma, sigma);	    
    
    coef = (v1 - v2) * v3;
  }
  return coef;
}

cx_mat calc_eff_Raman_operator(hoppings_bilayer2& ts, double delta, double kx, double ky, double kz, cx_double omega, double omega_i, cx_double *R){
  /* Mean field eigenenergies */
  cx_double ek1 = ts.ek1(kx, ky, kz);
  cx_double tz = ts.tz; 
  cx_double ek23 = ts.ek23(kx, ky, kz);
  cx_double ekz = ts.ekz(kx, ky, kz);	      
  double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, tz, kz, delta);
  double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, tz, kz, delta);

  /* Shifted omegas */
  // cx_double omega_i_shifted(omega_i, 0.5 * std::imag(omega));
  // cx_double omega_f_shifted(omega_i - std::real(omega), - 0.5 * std::imag(omega));
  double omega_f = omega_i - std::real(omega);
  
  /* Coefficients */
  // cx_double denom1 = ek_plus - ek_minus - omega_i_shifted;
  // cx_double denom2 = ek_plus - ek_minus + omega_f_shifted;  
  cx_double denom1 = ek_plus - ek_minus - omega_i;
  cx_double denom2 = ek_plus - ek_minus + omega_f;
  
  /* Formulation using up spin */
  cx_double x_k = xk(up_spin, ek1, tz, kz, delta);
  double z_k = zk(ek1, tz, kz, delta);
  
  cx_double uk = - 0.5 * x_k * z_k;
  cx_double vk = - 0.5 * x_k;
  cx_double wk = 0.5 * sqrt(1. - z_k * z_k);
  
  /* Effective Raman operator */
  mat tau_p = mat({0., 0., 1., 0.});
  mat tau_m = mat({0., 1., 0., 0.});
  mat tau_z = mat({1., 0., 0., - 1.});
  tau_p.set_size(2,2);
  tau_m.set_size(2,2);
  tau_z.set_size(2,2);

  cx_mat Pauli_0 = cx_mat({1., 0., 0., 1.});
  cx_mat Pauli_z = cx_mat({1., 0., 0., - 1.});
  Pauli_0.set_size(2,2);
  Pauli_z.set_size(2,2);

  cx_double rk = R[0] / denom1 + R[1] / denom2;
  cx_mat Mk = rk * (arma::kron(uk * tau_p + std::conj(uk) * tau_m, Pauli_0) + arma::kron(vk * tau_p - std::conj(vk) * tau_m + wk * tau_z, Pauli_z));
  return Mk;
}

cx_double calc_Raman_sigma0(int L, hoppings_bilayer2& ts, double ch_pot, double U, double T, double delta, CubaParam const& cbp, MatElemK const& me_K, cx_double omega, double omega_i, bool continuous_k){
  cx_double sum = 0;
  if ( continuous_k ) {
    std::cerr << "Function " << __func__ << " does not support continuous k.\n";
    std::exit(EXIT_FAILURE);
  } else { /* Integral for a finite-size system */
    double k1 = 2. * M_PI / L;
    for(int z=0; z < 2; ++z){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; ++x){    
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; ++y){
	  double ky = k1 * y;

	  /* Checking if the wavevector is inside the BZ. */
	  double factor = BZ_factor_square_half_filling(kx, ky);
	  if ( std::abs(factor) < 1e-12 ) { continue; }

	  /* Getting the matrix elements */
	  cx_double R[2];
	  if ( me_K.is_table_set() ) {
	    me_K.get_elem(kx, ky, kz, R);
	  } else {
	    me_K.calc_mat_elems(ts, delta, kx, ky, kz, R);
	  }

	  /* Effective Raman operator */	  
	  cx_mat Mk = calc_eff_Raman_operator(ts, delta, kx, ky, kz, omega, omega_i, R);	  

	  cx_double ek1 = ts.ek1(kx, ky, kz);
	  cx_double tz = ts.tz;	  
	  for(int sg1=-1; sg1<=1; sg1+=2){
	    for(int sg2=-1; sg2<=1; sg2+=2){
	      /* Prefactor */
	      cx_double prefactor = calc_prefactor_bare_res_func_bilayer(sg1, sg2, ts, T, kx, ky, kz, 0., 0., 0., omega, delta, ch_pot);
	      if ( std::abs(prefactor) < 1e-12 ) { continue; }
	      
	      /* Constructing projection operators */
	      cx_vec X1up = gs_HF1(up_spin, sg1, ek1, tz, kz, delta);
	      cx_vec X2up = gs_HF1(up_spin, sg2, ek1, tz, kz, delta);	      
	      cx_vec X1down = gs_HF1(down_spin, sg1, ek1, tz, kz, delta);
	      cx_vec X2down = gs_HF1(down_spin, sg2, ek1, tz, kz, delta);	      	      
	      cx_mat P1up = arma::kron(arma::trans(X1up), X1up);
	      cx_mat P2up = arma::kron(arma::trans(X2up), X2up);	      
	      cx_mat P1down = arma::kron(arma::trans(X1down), X1down);
	      cx_mat P2down = arma::kron(arma::trans(X2down), X2down);	      

	      /* Matrix element */
	      cx_double Kup = arma::trace(P1up * Mk * P2up * Mk.t());
	      cx_double Kdown = arma::trace(P1down * Mk * P2down * Mk.t());	      
	      cx_double K = Kup + Kdown;
	      	    
	      prefactor *= factor;
	      
	      sum += prefactor * K;
	    }
	  }
	}
      }
    }
    
    int n_units = L * L;  // Number of unit cells
    sum /= (double)(n_units);
  }
  
  return sum;
}

std::tuple<cx_double, cx_double> calc_Raman_sigma1(int L, hoppings_bilayer2& ts, double ch_pot, double U, double T, double delta, CubaParam const& cbp, MatElemK const& me_K, MatElemF const& me_F, cx_double omega, double omega_i, bool continuous_k){
  cx_double sigma1_00 = 0, sigma1_zz = 0;
  if ( continuous_k ) {
    std::cerr << "Function " << __func__ << " does not support continuous k.\n";
    std::exit(EXIT_FAILURE);
  } else { /* Integral for a finite-size system */
    double k1 = 2. * M_PI / L;

    /* Calculating the bubble of the internal vertices. */
    cx_mat chi0_00_up(NSUBL, NSUBL, arma::fill::zeros);
    cx_mat chi0_00_down(NSUBL, NSUBL, arma::fill::zeros);    
    cx_mat chi0_zz_up(NSUBL, NSUBL, arma::fill::zeros);
    cx_mat chi0_zz_down(NSUBL, NSUBL, arma::fill::zeros);        
    for(int z=0; z < 2; ++z){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; ++x){    
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; ++y){
	  double ky = k1 * y;

	  /* Checking if the wavevector is inside the BZ. */
	  double factor = BZ_factor_square_half_filling(kx, ky);
	  if ( std::abs(factor) < 1e-12 ) { continue; }
      
	  for(int sg1=-1; sg1<=1; sg1+=2){
	    for(int sg2=-1; sg2<=1; sg2+=2){
	      /* Prefactor */
	      cx_double prefactor = calc_prefactor_bare_res_func_bilayer(sg1, sg2, ts, T, kx, ky, kz, 0., 0., 0., omega, delta, ch_pot);
	      if ( std::abs(prefactor) < 1e-12 ) { continue; }
	      
	      int sg1i = (sg1+1) >> 1;    
	      int sg2i = (sg2+1) >> 1;

	      cx_double F00_up[NSUBL*NSUBL];
	      cx_double F00_down[NSUBL*NSUBL];      	      
	      cx_double Fzz_up[NSUBL*NSUBL];
	      cx_double Fzz_down[NSUBL*NSUBL];    	      

	      if ( me_F.is_table_set() ) {
		me_F.get_00_up(kx, ky, kz, sg1i, sg2i, F00_up);
		me_F.get_00_down(kx, ky, kz, sg1i, sg2i, F00_down);		
		me_F.get_zz_up(kx, ky, kz, sg1i, sg2i, Fzz_up);
		me_F.get_zz_down(kx, ky, kz, sg1i, sg2i, Fzz_down);		
	      } else {
		cx_double Fpm[NSUBL*NSUBL];		
		me_F.calc_mat_elems(ts, delta, kx, ky, kz, sg1, sg2, F00_up, F00_down, Fpm, Fzz_up, Fzz_down);
	      }

	      prefactor *= factor;

	      for(int i=0; i < 2; ++i){
		for(int j=0; j < 2; ++j){
		  int k = (i << 1) | j;
		  chi0_00_up(i,j) += prefactor * F00_up[k];
		  chi0_00_down(i,j) += prefactor * F00_down[k];
		  chi0_zz_up(i,j) += prefactor * Fzz_up[k];
		  chi0_zz_down(i,j) += prefactor * Fzz_down[k];
		}
	      }		  
	    }
	  }
	}
      }
    }
    int n_units = L * L;  // Number of unit cells
    chi0_00_up /= (double)(n_units);
    chi0_00_down /= (double)(n_units);
    chi0_zz_up /= (double)(n_units);
    chi0_zz_down /= (double)(n_units);
    cx_mat chi0_00 = chi0_00_up + chi0_00_down;    
    cx_mat chi0_zz = chi0_zz_up + chi0_zz_down;
    
    /* RPA */
    cx_mat denom_00 = arma::eye<arma::cx_mat>(NSUBL, NSUBL) + 0.5 * U * chi0_00;
    cx_mat denom_zz = arma::eye<arma::cx_mat>(NSUBL, NSUBL) - 0.5 * U * chi0_zz;    
    cx_mat U_eff_00 = - 0.5 * U * arma::inv(denom_00);   // Negative sign
    cx_mat U_eff_zz = 0.5 * U * arma::inv(denom_zz);
    
    /* Calculating the bubble of the internal and external vertices. */
    cx_vec Pi_0_up(NSUBL, arma::fill::zeros);
    cx_vec Pi_0_down(NSUBL, arma::fill::zeros);    
    cx_vec Pi_z_up(NSUBL, arma::fill::zeros);
    cx_vec Pi_z_down(NSUBL, arma::fill::zeros);            
    for(int z=0; z < 2; ++z){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; ++x){    
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; ++y){
	  double ky = k1 * y;

	  /* Checking if the wavevector is inside the BZ. */
	  double factor = BZ_factor_square_half_filling(kx, ky);
	  if ( std::abs(factor) < 1e-12 ) { continue; }

	  /* Getting the matrix elements */
	  cx_double Rup[2];
	  if ( me_K.is_table_set() ) {
	    me_K.get_elem(kx, ky, kz, Rup);
	  } else {
	    me_K.calc_mat_elems(ts, delta, kx, ky, kz, Rup);
	  }

	  /* Effective Raman operator */	  
	  cx_mat Mk = calc_eff_Raman_operator(ts, delta, kx, ky, kz, omega, omega_i, Rup);
	  
	  /* Sigma matrices */
	  mat tau_A = mat({1., 0., 0., 0});
	  mat tau_B = mat({0., 0., 0., 1});	  
	  tau_A.set_size(2,2);
	  tau_B.set_size(2,2);

	  cx_mat Pauli_0 = cx_mat({1., 0., 0., 1.});
	  cx_mat Pauli_z = cx_mat({1., 0., 0., - 1.});
	  Pauli_0.set_size(2,2);
	  Pauli_z.set_size(2,2);

	  cx_mat Sigma_A0 = arma::kron(tau_A, Pauli_0);
	  cx_mat Sigma_B0 = arma::kron(tau_B, Pauli_0);
	  cx_mat Sigma_Az = arma::kron(tau_A, Pauli_z);
	  cx_mat Sigma_Bz = arma::kron(tau_B, Pauli_z);
	  
	  for(int sg1=-1; sg1<=1; sg1+=2){
	    for(int sg2=-1; sg2<=1; sg2+=2){
	      /* Prefactor */
	      cx_double prefactor = calc_prefactor_bare_res_func_bilayer(sg1, sg2, ts, T, kx, ky, kz, 0., 0., 0., omega, delta, ch_pot);
	      if ( std::abs(prefactor) < 1e-12 ) { continue; }
	      
	      /* Constructing projection operators */
	      cx_double ek1 = ts.ek1(kx, ky, kz);
	      cx_double tz = ts.tz;	      
	      cx_vec X1up = gs_HF1(up_spin, sg1, ek1, tz, kz, delta);
	      cx_vec X2up = gs_HF1(up_spin, sg2, ek1, tz, kz, delta);
	      cx_vec X1down = gs_HF1(down_spin, sg1, ek1, tz, kz, delta);
	      cx_vec X2down = gs_HF1(down_spin, sg2, ek1, tz, kz, delta);	      
	      cx_mat P1up = arma::kron(arma::trans(X1up), X1up);
	      cx_mat P2up = arma::kron(arma::trans(X2up), X2up);	      
	      cx_mat P1down = arma::kron(arma::trans(X1down), X1down);
	      cx_mat P2down = arma::kron(arma::trans(X2down), X2down);	      

	      /* Matrix element */
	      cx_double P_0_A_up = arma::trace(P1up * Sigma_A0 * P2up * Mk.t());
	      cx_double P_0_B_up = arma::trace(P1up * Sigma_B0 * P2up * Mk.t());
	      cx_double P_0_A_down = arma::trace(P1down * Sigma_A0 * P2down * Mk.t());
	      cx_double P_0_B_down = arma::trace(P1down * Sigma_B0 * P2down * Mk.t());	      
	      cx_double P_z_A_up = arma::trace(P1up * Sigma_Az * P2up * Mk.t());
	      cx_double P_z_B_up = arma::trace(P1up * Sigma_Bz * P2up * Mk.t());
	      cx_double P_z_A_down = arma::trace(P1down * Sigma_Az * P2down * Mk.t());
	      cx_double P_z_B_down = arma::trace(P1down * Sigma_Bz * P2down * Mk.t());	      
	      
	      prefactor *= factor;

	      Pi_0_up(0) += prefactor * P_0_A_up;
	      Pi_0_up(1) += prefactor * P_0_B_up;
	      Pi_0_down(0) += prefactor * P_0_A_down;
	      Pi_0_down(1) += prefactor * P_0_B_down;	      
	      Pi_z_up(0) += prefactor * P_z_A_up;
	      Pi_z_up(1) += prefactor * P_z_B_up;
	      Pi_z_down(0) += prefactor * P_z_A_down;
	      Pi_z_down(1) += prefactor * P_z_B_down;
	    }
	  }
	}
      }
    }    
    Pi_0_up /= (double)(n_units);
    Pi_0_down /= (double)(n_units);    
    Pi_z_up /= (double)(n_units);
    Pi_z_down /= (double)(n_units);
    cx_vec Pi_0 = Pi_0_up + Pi_0_down;    
    cx_vec Pi_z = Pi_z_up + Pi_z_down;
    
    /* Total contributions */
    cx_mat sigma1_00_ = Pi_0.st() * U_eff_00 * arma::conj(Pi_0);
    sigma1_00 = sigma1_00_[0];
    cx_mat sigma1_zz_ = Pi_z.st() * U_eff_zz * arma::conj(Pi_z);
    sigma1_zz = sigma1_zz_[0];
  }
  
  return std::make_tuple(sigma1_00, sigma1_zz);
}

void calc_Raman_bilayer(path& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  int Lk = pr.Lk;
  double U = pr.U;
  double filling = pr.filling;
  double T = pr.T;
  bool continuous_k = pr.continuous_k;

  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Omegas */
  double omega_min = pr.omega_min;
  double omega_max = pr.omega_max;
  double omega_delta = pr.omega_delta;  
  
  int n_omegas = int((omega_max - omega_min)/omega_delta+0.5);
  std::vector<double> omegas(n_omegas);
  for(int o=1; o <= n_omegas; ++o){ omegas[o-1] = omega_min + omega_delta * o; }
  
  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  /* Calculate the chemical potential and the charge gap. */
  double delta = solve_self_consistent_eq_bilayer2( L, *ts, U, filling, T, cbp, continuous_k );  
  std::cout << "delta = " << delta << std::endl;
  /* Assume that the chemical potential does not depend on L for integral over continuous k. */
  double ch_gap, ch_pot;
  std::tie(ch_gap, ch_pot) = calc_charge_gap_bilayer( L, *ts, delta );  /* Finite size */  

  /* Output */
  ofstream out_gap;
  out_gap.open(base_dir / "gap.text");
  out_gap << "Charge gap = " << ch_gap << std::endl;
  out_gap << "Chemical potential = " << ch_pot << std::endl;
  
  /* MatElemF */
  MatElemF me_F( L, L, 2, NSUBL );

  auto calc_gaps = [&](double qx, double qy, double qz){
    /* Setting the ordering vector */
    me_F.set_q(qx, qy, qz);
    if ( !continuous_k ) {
      /* The polarizations are calculated in advance. */
      me_F.set_table( *ts, delta );
    }
    
    /* Calculating gaps */
    double omega_T = 0, omega_L = 0, omega_ph = 0;
    bool return_upper = true;
    bool verbose = false;    
    std::tie(omega_T, omega_L, omega_ph) = calc_gap_bilayer(L, *ts, ch_pot, U, T, delta, cbp, me_F, continuous_k, return_upper, verbose);
    
    return std::make_tuple(omega_T, omega_L, omega_ph);
  };

  /* Calculating gaps */
  double omega_T = 0, omega_L = 0, omega_ph = 0;
  std::tie(omega_T, omega_L, omega_ph) = calc_gaps(0., 0., 0.);
  
  /* Output */
  out_gap << "omega_T = " << omega_T << std::endl;
  out_gap << "omega_L = " << omega_L << std::endl;
  out_gap << "omega_ph = " << omega_ph << std::endl;
  
  /* Calculating Psi'*/
  // double Psider = 1.0;
  // std::cout << "Setting Psider to one for simplicity." << std::endl;
  double Psider = calc_Psider(L, *ts, ch_pot, omega_L, delta, cbp, continuous_k);

  /* Output */
  out_gap << "Psider = " << Psider << std::endl;

  /* Calculating the minimum bk square */
  double min_bk_sq = calc_min_bk_sq_bilayer(L, *ts);
    
  /* Constants */
  // double stag_rot_angle = 12.0 / 180.0 * M_PI;
  double k1 = 2. * M_PI / (double)L;  
  double g_photon = 1.0;  // This factor does not matter to the result.
  // double g_photon = sqrt(planck_h * c_light * c_light / (pr.omega_i * 2. * M_PI / planck_h));  
  double beta = 0;
  if ( !pr.T_equal_to_0 ) {
    beta = 1. / (kB * pr.T);
  }

  /* Precision */
  int prec = 10;
  int pw = prec + 10;

  // /* Calculating the scattering states of a single particle-hole pair. */
  // int N = 2 * L * L;
  // vec vals(N);
  // cx_mat Uph1(N, N);
  // calc_particle_hole1(pr, *ts, delta, vals, Uph1);
      
  std::vector<BondDelta> bonds {
     {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, -1, 0},  {0, 0, 1}     
   };

  // auto wavefunc_sublattice_correction = [](int g1, int g2, cx_double psi){
  //   if (g1 == 0) {
  //     if (g2 == 0) {
  // 	return psi;
  //     } else {
  // 	return psi;
  //     }
  //   } else {
  //     if (g2 == 0) {
  // 	return std::conj(psi);
  //     } else {
  // 	return - psi;
  //     }
  //   }
  // };
  
  // /* Output of the mean field energy gaps. */
  // ofstream out_MF_gap;
  // out_MF_gap.open(base_dir / "mean_field_energy_gap.text");    
  // out_MF_gap << "# Index  kx  ky  kz  E_plus  E_minus  E_gap" << std::endl;  

  // std::vector<double> particle_hole_gaps;
  // for(int z=0; z < 2; ++z){    
  //   double kz = M_PI * z;
  //   for(int x=-L/2; x < L/2; ++x){
  //     double kx = k1 * x;
  //     for(int y=-L/2; y < L/2; ++y){
  // 	double ky = k1 * y;
	  
  // 	/* Checking if the wavevector is inside the BZ. */
  // 	double factor = BZ_factor_square_half_filling(kx, ky);
  // 	if ( std::abs(factor) < 1e-12 ) { continue; }

  // 	/* Eigenenergy */
  // 	cx_double ek1 = ts->ek1(kx, ky, kz);
  // 	cx_double ek23 = ts->ek23(kx, ky, kz);
  // 	cx_double ekz = ts->ekz(kx, ky, kz);	      
  // 	double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts->tz, kz, delta);
  // 	double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts->tz, kz, delta);
  // 	particle_hole_gaps.push_back(ek_plus - ek_minus);
	
  // 	out_MF_gap << kx << "  " << ky << "  " << kz << "  " << ek_plus << "  " << ek_minus << "  " << ek_plus - ek_minus << std::endl;
  //     }
  //   }
  // }

  // /* Sorted mean-field particle-hole energy gaps */
  // std::sort(particle_hole_gaps.begin(), particle_hole_gaps.end());
  // ofstream out_eh_gap(base_dir / "particle_hole_gaps.text");
  // double e = 0;
  // for(double eg: particle_hole_gaps){
  //   if ( std::abs(eg - e) > 1e-10 ) {
  //     out_eh_gap << eg << std::endl;
  //     e = eg;
  //   }
  // }
  // out_eh_gap.close();
  
  // /* Output of the dispersion calculated by RPA. */  
  // ofstream out_dispersion_transverse, out_dispersion_longitudinal;
  // out_dispersion_transverse.open(base_dir / "dispersion_transverse.text");
  // out_dispersion_longitudinal.open(base_dir / "dispersion_longitudinal.text");
  // out_dispersion_transverse << "# Index  kx  ky  kz  gap" << std::endl;
  // out_dispersion_longitudinal << "# Index  kx  ky  kz  gap" << std::endl;
  
  // /* Output of the wave function. */
  // ofstream out_wavefunc_k;
  // out_wavefunc_k.open(base_dir / "wavefunction_k.text");
  // out_wavefunc_k << "# Index  kx  ky  kz   Re[psi_AA]   Im[psi_AA]   Re[psi_AB]   Im[psi_AB]" << std::endl;
  
  // /* Energies. */
  // std::vector<double> E_k_plus;
  // std::vector<double> E_k_minus;
  // std::vector<double> gap_T;
  // std::vector<double> gap_L;

  // /* Setting the parameters */
  // int spin = up_spin;  // Unused  
  // int sublattice = 0;  // Unused
  // int diff_r[] = {0,0,0};  // Unused

  // /* Photon momenta */
  // std::vector<vec3> photon_Qs;
  // photon_Qs.push_back(vec3{0, 0, 0});    
  // photon_Qs.push_back(vec3{   M_PI,   M_PI,   M_PI });
  // // photon_Qs.push_back(vec3{ - M_PI, - M_PI, - M_PI });
  // // std::vector<std::pair<int,int>> photon_ks{{1,0}};
  // std::vector<std::pair<int,int>> photon_ks{{0,0}};
  // // std::vector<std::pair<int,int>> photon_ks{{1,1}};    
  // // std::vector<std::pair<int,int>> photon_ks{{0,0},{0,1},{1,0},{1,1}};

  // std::size_t idx_k_sigma_ki_kf = 0;
  // for(int z=0; z < 2; ++z){    
  //   double kz = M_PI * z;
  //   for(int x=-L/2; x < L/2; ++x){
  //     double kx = k1 * x;
  //     for(int y=-L/2; y < L/2; ++y){
  // 	double ky = k1 * y;
	  
  // 	/* Checking if the wavevector is inside the BZ. */
  // 	double factor = BZ_factor_square_half_filling(kx, ky);
  // 	if ( std::abs(factor) < 1e-12 ) { continue; }

  // 	/* Calculating the gaps */
  // 	double gT = 0, gL = 0, gph = 0;
  // 	std::tie(gT, gL, gph) = calc_gaps(kx, ky, kz);  // Solving the pole equation in RPA.
	    
  // 	for(int sigma: {up_spin, down_spin}){
  // 	  for(std::pair<int,int> photon_k: photon_ks){	  
  // 	    gap_T.push_back(gT);
  // 	    gap_L.push_back(gL);

  // 	    /* Eigenenergy */
  // 	    cx_double ek1 = ts->ek1(kx, ky, kz);
  // 	    cx_double ek23 = ts->ek23(kx, ky, kz);
  // 	    cx_double ekz = ts->ekz(kx, ky, kz);	      
  // 	    double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts->tz, kz, delta);
  // 	    double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts->tz, kz, delta);
  // 	    E_k_plus.push_back(ek_plus);
  // 	    E_k_minus.push_back(ek_minus);

  // 	    /* For checking the wave functions. */
  // 	    wfib.set_parameters(*ts, pr.largeUlimit, pr.largeU_scaling_prefactor, sigma, omega_L, Psider, delta, diff_r, sublattice, min_bk_sq);	  
  // 	    cx_double xki = xk(wfib.spin(), ek1, ts->tz, kz, wfib.delta());
  // 	    double zki = zk(ek1, ts->tz, kz, wfib.delta());
  // 	    cx_double bki = bk(wfib.spin(), ek1, ts->tz, kz);
  // 	    cx_double wavefunc_AA = wfib.integrand(xki, zki, bki, ek_plus, ek_minus, 1., 0);
  // 	    cx_double wavefunc_AB = wfib.integrand(xki, zki, bki, ek_plus, ek_minus, 1., 1);	
		    
  // 	    /* Output */
  // 	    int npw = 12;
  // 	    out_dispersion_transverse << idx_k_sigma_ki_kf << std::setw(npw) << kx << std::setw(npw) << ky << std::setw(npw) << kz << std::setw(pw) << gT << std::endl;
  // 	    out_dispersion_longitudinal << idx_k_sigma_ki_kf << std::setw(npw) << kx << std::setw(npw) << ky << std::setw(npw) << kz << std::setw(pw) << gL << std::endl;

  // 	    out_wavefunc_k << idx_k_sigma_ki_kf << std::setw(npw) << kx << std::setw(npw) << ky << std::setw(npw) << kz << std::setw(pw) << std::real(wavefunc_AA) << std::setw(pw) << std::imag(wavefunc_AA) << std::setw(pw) << std::real(wavefunc_AB) << std::setw(pw) << std::imag(wavefunc_AB) << std::endl;;
	
  // 	    ++idx_k_sigma_ki_kf;
  // 	  }
  // 	}
  //     }
  //   }
  // }

  auto invalid_components = [](int mu, int nu){
    /* Only considering in-plane components */
    if (mu == 2 || nu == 2) { return true; }
    else { return false; }
    
    // /* Not considering xz, yz, zx, and zy. */
    // if (mu == 2 && nu != 2) { return true; }
    // else if ( mu != 2 && nu == 2) { return true; }
    // else { return false; }
  };
  
  /* Results */
  int n_spec = 3;   // Number of spectrum types
  mat spec_Raman_xx(n_spec, n_omegas);
  mat spec_Raman_xy(n_spec, n_omegas);
  mat spec_Raman_xz(n_spec, n_omegas);
  mat spec_Raman_yx(n_spec, n_omegas);
  mat spec_Raman_yy(n_spec, n_omegas);
  mat spec_Raman_yz(n_spec, n_omegas);
  mat spec_Raman_zx(n_spec, n_omegas);
  mat spec_Raman_zy(n_spec, n_omegas);
  mat spec_Raman_zz(n_spec, n_omegas);    

  auto spec_Raman = [&](int mu, int nu){
    mat *spec;
    if ( mu == 0 ) {
      if ( nu == 0 ) { spec = &spec_Raman_xx; }
      else if ( nu == 1 ) { spec = &spec_Raman_xy; }
      else { spec = &spec_Raman_xz; }
    } else if ( mu == 1 ) {
      if ( nu == 0 ) { spec = &spec_Raman_yx; }
      else if ( nu == 1 ) { spec = &spec_Raman_yy; }
      else { spec = &spec_Raman_yz; }
    } else {
      if ( nu == 0 ) { spec = &spec_Raman_zx; }
      else if ( nu == 1 ) { spec = &spec_Raman_zy; }
      else { spec = &spec_Raman_zz; }      
    }
    return spec;
  };
	      
  // /* Calculating the spectral weights and the energies. */    
  // // std::vector<std::vector<cx_double>> weight_plus;
  // std::vector<std::vector<cx_double>> weight_minus1, weight_minus2;
  // for(int mu=0; mu < 3; ++mu){
  //   BondDelta e_mu(mu);
  //   for(int nu=0; nu < 3; ++nu){
  //     BondDelta e_nu(nu);
      
  //     // std::vector<cx_double> weight_plus_comp;
  //     std::vector<cx_double> weight_minus1_comp, weight_minus2_comp;      
  //     idx_k_sigma_ki_kf = 0;
  //     for(int z=0; z < 2; ++z){    
  // 	double kz = M_PI * z;
  // 	for(int x=-L/2; x < L/2; ++x){
  // 	  double kx = k1 * x;
  // 	  for(int y=-L/2; y < L/2; ++y){
  // 	    double ky = k1 * y;

  // 	    /* Checking if the wavevector is inside the BZ. */
  // 	    double factor = BZ_factor_square_half_filling(kx, ky);
  // 	    if ( std::abs(factor) < 1e-12 ) { continue; }
	    
  // 	    cx_double ek1 = ts->ek1(kx, ky, kz);
  // 	    cx_mat Uk = gs_HF(ek1, ts->tz, kz, delta);
  // 	    cx_mat Uk_bar = gs_HF(ek1, ts->tz, kz + M_PI, delta);	    
  // 	    cx_mat Uk_dg = Uk.t();

  // 	    for(int sigma: {up_spin, down_spin}){  // Electron spin in the conduction band labeling the final state
  // 	      for(std::pair<int,int> photon_k: photon_ks){
  // 		/* Velocity */
  // 		int kii = photon_k.first;
  // 		int kfi = photon_k.second;
		
  // 		vec3 ki = photon_Qs[kii];
  // 		vec3 kf = photon_Qs[kfi];
		    
  // 		cx_double v_b2_34_pp = velocity_U1(*ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, kf, e_nu, bonds, 1, 1, sigma, sigma);
  // 		cx_double v_b1_12_pm = velocity_U1(*ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, - ki, e_mu, bonds, 1, -1, sigma, sigma);
  // 		cx_double v_b1_12_pp = velocity_U1(*ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, - ki, e_mu, bonds, 1,  1, sigma, sigma);
  // 		cx_double v_b2_34_pm = velocity_U1(*ts, Uk_dg, Uk, Uk_bar, kx, ky, kz, kf, e_nu, bonds, 1, -1, sigma, sigma);
			
  // 		cx_double v1 = v_b2_34_pp * v_b1_12_pm;
  // 		cx_double v2 = v_b1_12_pp * v_b2_34_pm;
		
  // 		/* Matrix elements */
  // 		cx_double coef = 0.0;
  // 		for(int g5=0; g5 < 2; ++g5){
  // 		  int g5_idx = sublattice_spin_index(g5, sigma);
  // 		  cx_double U_5plus = Uk(g5_idx, sign_spin_index(1, sigma));
		  
  // 		  for(int g6=0; g6 < 2; ++g6){
  // 		    int g6_idx = sublattice_spin_index(g6, sigma);
  // 		    // cx_double U6_plus = Uk(g6_idx, sign_spin_index(1, sigma));
  // 		    cx_double Udg_minus6 = Uk_dg(sign_spin_index(-1, sigma), g6_idx);
		    
  // 		    /* Using the flipped spin for the wave-function calculation depending on the sublattices. */
  // 		    int sigma_psi = sigma;
  // 		    if ( g6 == 1 && g5 == 0 ) {
  // 		      if ( sigma == up_spin ) {
  // 			sigma_psi = down_spin;
  // 		      } else {
  // 			sigma_psi = up_spin;
  // 		      }
  // 		    }

  // 		    /* Setting the parameters. */
  // 		    int diff_r[] = {0,0,0};  // Unused
  // 		    int sublattice = g5 == g6 ? 0 : 1;		    
  // 		    wfib.set_parameters(*ts, pr.largeUlimit, pr.largeU_scaling_prefactor, sigma_psi, omega_L, Psider, delta, diff_r, sublattice, min_bk_sq);
		
  // 		    /* Calculating Psi(g6, sigma, g5, sigma) from Psi(A, sigma_psi, A/B, sigma_psi) */
  // 		    cx_double xki = xk(wfib.spin(), ek1, ts->tz, kz, wfib.delta());
  // 		    double zki = zk(ek1, ts->tz, kz, wfib.delta());
  // 		    cx_double bki = bk(wfib.spin(), ek1, ts->tz, kz);

  // 		    /* Mean field eigenenergy */
  // 		    double ek_plus = E_k_plus[idx_k_sigma_ki_kf];
  // 		    double ek_minus = E_k_minus[idx_k_sigma_ki_kf];
		  
  // 		    /* Psi(A, sigma_psi, A/B, sigma_psi) */
  // 		    cx_double wavefunc = wfib.integrand(xki, zki, bki, ek_plus, ek_minus, 1., sublattice);

  // 		    /* Transforming to Psi(g6, sigma, g5, sigma) */
  // 		    wavefunc = wavefunc_sublattice_correction(g6, g5, wavefunc);

  // 		    /* Adding the matrix elements. */
  // 		    coef += U_5plus * Udg_minus6 * std::conj(wavefunc);
  // 		  }  /* end for g6 */	    
  // 		}  /* end for g5 */
	      
  // 		cx_double mat_elem1_minus = coef * v1;
  // 		cx_double mat_elem2_minus = coef * v2;
		  
  // 		weight_minus1_comp.push_back(mat_elem1_minus);
  // 		weight_minus2_comp.push_back(mat_elem2_minus);

  // 		++idx_k_sigma_ki_kf;
  // 	      } /* end for photon_k */	      
  // 	    } /* end for sigma */	      
  // 	  } /* end for y */
  // 	} /* end for x */
  //     } /* end for z */
  
  // 	// weight_plus.push_back(weight_plus_comp);
  //     weight_minus1.push_back(weight_minus1_comp);
  //     weight_minus2.push_back(weight_minus2_comp);      
  //   } /* end for nu */
  // } /* end for mu */

  /* Matrix elements for the susceptibility */
  me_F.set_q(0., 0., 0.);
  me_F.set_table(*ts, delta);

  // // for check
  // cx_double uk(1., 2.);
  // cx_double vk(10., 20);
  // cx_double wk(100., 0.);
  
  // mat tau_p = mat({0., 0., 1., 0.});
  // mat tau_m = mat({0., 1., 0., 0.});
  // mat tau_z = mat({1., 0., 0., - 1.});
  // tau_p.set_size(2,2);
  // tau_m.set_size(2,2);
  // tau_z.set_size(2,2);

  // cx_mat Pauli_0 = cx_mat({1., 0., 0., 1.});
  // cx_mat Pauli_z = cx_mat({1., 0., 0., - 1.});
  // Pauli_0.set_size(2,2);
  // Pauli_z.set_size(2,2);

  // cx_double rk = 10.0;
  // cx_mat Mk = rk * (arma::kron(uk * tau_p + std::conj(uk) * tau_m, Pauli_0) + arma::kron(vk * tau_p - std::conj(vk) * tau_m + wk * tau_z, Pauli_z));
  // // cx_mat Mk = arma::kron(uk * tau_p + std::conj(uk) * tau_m, Pauli_0) + arma::kron(vk * tau_p - std::conj(vk) * tau_m + wk * tau_z, Pauli_z);
  
  // std::cout << "Test Mk:" << std::endl;
  // std::cout << Mk << std::endl;
  // std::cout << Mk.t() << std::endl;
  // std::cout << Mk << std::endl;

  
  /* Calculating the Raman scattering cross sections. */
#ifdef WITH_OpenMP
#pragma omp parallel
  {
  /* Assignment for each thread */  
  int n_threads = omp_get_num_threads();
  int nt = n_omegas / n_threads;
  int rem = n_omegas % n_threads;      
  int thread_id = omp_get_thread_num();
  if ( thread_id < rem ) { nt += 1; }
      
  /* Results for each thread */
  mat spec_Raman_xx_thread(n_spec, nt);
  mat spec_Raman_xy_thread(n_spec, nt);
  mat spec_Raman_xz_thread(n_spec, nt);
  mat spec_Raman_yx_thread(n_spec, nt);
  mat spec_Raman_yy_thread(n_spec, nt);
  mat spec_Raman_yz_thread(n_spec, nt);
  mat spec_Raman_zx_thread(n_spec, nt);
  mat spec_Raman_zy_thread(n_spec, nt);
  mat spec_Raman_zz_thread(n_spec, nt);      

  auto spec_Raman_thread = [&](int mu, int nu){
    mat *spec;
    if ( mu == 0 ) {
      if ( nu == 0 ) { spec = &spec_Raman_xx_thread; }
      else if ( nu == 1 ) { spec = &spec_Raman_xy_thread; }
      else { spec = &spec_Raman_xz_thread; }
    } else if ( mu == 1 ) {
      if ( nu == 0 ) { spec = &spec_Raman_yx_thread; }
      else if ( nu == 1 ) { spec = &spec_Raman_yy_thread; }
      else { spec = &spec_Raman_yz_thread; }
    } else {
      if ( nu == 0 ) { spec = &spec_Raman_zx_thread; }
      else if ( nu == 1 ) { spec = &spec_Raman_zy_thread; }
      else { spec = &spec_Raman_zz_thread; }      
    }
    return spec;
  };
  
  for(int oidx=0; oidx < nt; ++oidx){
    int o = thread_id + oidx * n_threads; // stride: n_threads
#else
  for(int o=0; o < n_omegas; ++o){
#endif    
    /* Energy of the final photon state */
    // double omega_f = pr.omega_i - omegas[o]; // (eV)
    cx_double omega_shifted = cx_double(omegas[o], pr.eta);
    for(int mu=0; mu < 3; ++mu){
      for(int nu=0; nu < 3; ++nu){
	if (invalid_components(mu,nu)) { continue; }	
	/* Matrix elements */
	MatElemK me_K(L, L, 2, mu, nu, bonds);

	/* Initialization */
	me_K.set_q(0., 0., 0.);
	me_K.set_table(*ts, delta);
	  
	cx_double sigma0 = calc_Raman_sigma0(L, *ts, ch_pot, U, T, delta, cbp, me_K, omega_shifted, pr.omega_i, continuous_k);
	cx_double sigma1_00, sigma1_zz;
	std::tie(sigma1_00, sigma1_zz) = calc_Raman_sigma1(L, *ts, ch_pot, U, T, delta, cbp, me_K, me_F, omega_shifted, pr.omega_i, continuous_k);

#ifdef WITH_OpenMP
	(*spec_Raman_thread(mu,nu))(0, oidx) = 2.0 * std::imag(sigma0);
	(*spec_Raman_thread(mu,nu))(1, oidx) = 2.0 * std::imag(sigma1_00);
	(*spec_Raman_thread(mu,nu))(2, oidx) = 2.0 * std::imag(sigma1_zz);	
#else
	(*spec_Raman(mu,nu))(0, o) = 2.0 * std::imag(sigma0);
	(*spec_Raman(mu,nu))(1, o) = 2.0 * std::imag(sigma1_00);
	(*spec_Raman(mu,nu))(2, o) = 2.0 * std::imag(sigma1_zz);		
#endif
      }
    }
  }
      
#ifdef WITH_OpenMP
  /* Extracting the results */
  for(int mu=0; mu < 3; ++mu){
    for(int nu=0; nu < 3; ++nu){
      if (invalid_components(mu,nu)) { continue; }
      for(int oidx=0; oidx < nt; ++oidx){
	int o = thread_id + oidx * n_threads;
	for(int ispec=0; ispec < n_spec; ++ispec){
	  (*spec_Raman(mu,nu))(ispec, o) = (*spec_Raman_thread(mu,nu))(ispec, oidx);
	}
      }
    }
  }
  }
#endif      

	// arma::cx_mat chi0_00(NSUBL,NSUBL,arma::fill::zeros);
	// arma::cx_mat chi0_zz(NSUBL,NSUBL,arma::fill::zeros);

	// /* Calculating the internal vertex. */
	// std::tie(chi0_00, chi0_zz) = calc_internal_bubble_bilayer(L, ts, mu, U, T, delta, cbp, me_F, omega_shifted, continuous_k);

	// arma::cx_mat Pi_0(NSUBL,1,arma::fill::zeros);
	// arma::cx_mat Pi_z(NSUBL,1,arma::fill::zeros);

	// /* Calculating the external vertex. */
	// std::tie(Pi_0, Pi_z) = calc_external_vertex_bilayer(L, ts, mu, U, T, delta, cbp, me_F, omega_shifted, continuous_k);
      
	// /* RPA */
	// /* Transverse = < \sigma^- \sigma^+ >; Longitudinal (zz) = < \sigma^z \sigma^z > */
	// /* Note that 2 < \sigma^- \sigma^+ > = < \sigma^z \sigma^z > (U -> 0 for the SU(2) case) */
	// arma::cx_mat denom_00 = arma::eye<arma::cx_mat>(NSUBL,NSUBL) + 0.5 * U * chi0_00;        
	// arma::cx_mat denom_zz = arma::eye<arma::cx_mat>(NSUBL,NSUBL) - 0.5 * U * chi0_zz;  
	// arma::cx_mat chi_00 = arma::inv(denom_pm) * chi0_00;
	// arma::cx_mat chi_zz = arma::inv(denom_zz) * chi0_zz;
	
	// // sigma-to-spin factor
	// cx_double chi_xy = 0.5 * arma::accu(chi_pm);
	// cx_double chi_z = 0.25 * arma::accu(chi_zz);

	// cx_double sigma_omega0 = calc_Raman_sigma0();
	// cx_double sigma_omega1 = calc_Raman_sigma1();	
	// cx_double sigma_omega = sigma_omega0 + sigma_omega1;


	
	// idx_k_sigma_ki_kf = 0;
	// cx_double sum = 0.0;	
	// for(int z=0; z < 2; ++z){    
	//   double kz = M_PI * z;
	//   for(int x=-L/2; x < L/2; ++x){
	//     double kx = k1 * x;
	//     for(int y=-L/2; y < L/2; ++y){
	//       double ky = k1 * y;

	//       /* Checking if the wavevector is inside the BZ. */
	//       double factor = BZ_factor_square_half_filling(kx, ky);
	//       if ( std::abs(factor) < 1e-12 ) { continue; }

	//       double e_gap = gap_L[idx_k_sigma_ki_kf];
	      
	//       cx_double weight_minus_sum = 0;
	//       double down_sign = 1.0;   // Linear combination

	//       /* Taking the sum over spins and photon momenta. */
	//       for(int sigma: {up_spin, down_spin}){
	// 	for(std::pair<int,int> photon_k: photon_ks){
	// 	  /* Factor */
	// 	  int kii = photon_k.first;
	// 	  int kfi = photon_k.second;		  
	// 	  double factor_q = kii - kfi == 0 ? 1.0 : 0.5;
		  
	// 	  /* Mean field eigenenergy */
	// 	  double ek_plus = E_k_plus[idx_k_sigma_ki_kf];
	// 	  double ek_minus = E_k_minus[idx_k_sigma_ki_kf];		  
	// 	  double ek_diff = ek_plus - ek_minus;
		
	// 	  cx_double weight1_minus = weight_minus1[idx_comp][idx_k_sigma_ki_kf] / (ek_diff - pr.omega_i);	      
	// 	  cx_double weight2_minus = weight_minus2[idx_comp][idx_k_sigma_ki_kf] / (ek_diff + omega_f);
	// 	  weight_minus_sum += factor_q * (weight1_minus + weight2_minus) * (sigma == up_spin ? 1.0 : down_sign) / sqrt(2.0);
	// 	  ++idx_k_sigma_ki_kf;
	// 	}  /* end for photon_k */
	//       }  /* end for sigma */

	//       // Zero temperature
	//       double Ei = 0;
	//       sum += factor * exp(- beta * Ei) * std::norm(g_photon * weight_minus_sum) / (e_gap - omega_shifted);	      
	//     } /* end for y */
	//   } /* end for x */
	// } /* end for z */
	
	// spec_Raman(mu,nu)[o] = 2.0 * std::imag(sum);
  
  // 	++idx_comp;
  //     } /* end for nu */
  //   } /* end for mu */    
  // }/* end for o */

  /* Output */
  for(int mu=0; mu < 3; ++mu){
    std::string mu_str;
    if ( mu == 0 ) { mu_str = 'x'; }
    else if ( mu == 1 ) { mu_str = 'y'; }
    else { mu_str = 'z'; }
    
    for(int nu=0; nu < 3; ++nu){
      if (invalid_components(mu,nu)) { continue; }      
      std::string nu_str;
      if ( nu == 0 ) { nu_str = 'x'; }
      else if ( nu == 1 ) { nu_str = 'y'; }
      else { nu_str = 'z'; }
  
      ofstream out_raman;
      std::string ofilen("Raman_scattering-"+mu_str+nu_str+".text");
      out_raman.open(base_dir / ofilen);      
      out_raman << "# Omega     sigma0     sigma1_00     sigma1_zz     sum" << std::endl;
  
      /* Output */
      mat *spec = spec_Raman(mu,nu);
      for(int o=0; o < n_omegas; ++o){
	double omega = omegas[o];
	out_raman << omega;
	
	double sum = 0;
	for(int ispec=0; ispec < n_spec; ++ispec){
	  double val = (*spec)(ispec, o);
	  sum += val;
	  out_raman << std::setw(pw) << val;
	}
	out_raman << std::setw(pw) << sum << std::endl;	
      } /* end for o */
  
      out_raman.close();        
    } /* end for nu */
  } /* end for mu */
  
  // /* Closing the output files. */
  // out_gap.close();
  // out_dispersion_transverse.close();
  // out_dispersion_longitudinal.close();
  // out_MF_gap.close();
  // out_wavefunc_k.close();
}

void calc_coef_eff_Raman_real_space(path& base_dir, rpa::parameters const& pr){
  /* Output of the coefficient of the effective Raman operator. */

  /* Getting parameters */
  int L = pr.L;
  double U = pr.U;
  double filling = pr.filling;
  double T = pr.T;
  bool continuous_k = pr.continuous_k;

  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);
  
  /* Parameters for Cuba */
  CubaParam cbp(pr);

  /* Calculate the chemical potential and the charge gap. */
  double delta = solve_self_consistent_eq_bilayer2(L, *ts, U, filling, T, cbp, continuous_k);  
  std::cout << "delta = " << delta << std::endl;
  
  /* Bonds */
  std::vector<BondDelta> bonds {
     {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, -1, 0},  {0, 0, 1}     
   };
  
  /* Shifted omegas */
  cx_double omega_shifted = cx_double(pr.Omega, pr.eta);
  cx_double omega_i_shifted(pr.omega_i, 0.5 * std::imag(omega_shifted));
  cx_double omega_f_shifted(pr.omega_i - std::real(omega_shifted), - 0.5 * std::imag(omega_shifted));
  // double omega_f = omega_i - std::real(omega);

  
  /* For each component set */
  for(int mu=0; mu < 3; ++mu){
    std::string mu_str;
    if ( mu == 0 ) { mu_str = 'x'; }
    else if ( mu == 1 ) { mu_str = 'y'; }
    else { mu_str = 'z'; }
    
    for(int nu=0; nu < 3; ++nu){
      std::string nu_str;
      if ( nu == 0 ) { nu_str = 'x'; }
      else if ( nu == 1 ) { nu_str = 'y'; }
      else { nu_str = 'z'; }

      ofstream out_coef_Raman;
      std::string ofilen("coefficient_of_Raman_operator-"+mu_str+nu_str+".text");
      out_coef_Raman.open(base_dir/ofilen);      
      out_coef_Raman << "# Coefficient of the Sz-like term for Omega = " << pr.Omega << " and omega_i = " << pr.omega_i << std::endl;      
      out_coef_Raman << "#  x   y   z     Re     Im" << std::endl;
	
      /* Matrix elements */
      MatElemK me_K(L, L, 2, mu, nu, bonds);
      me_K.set_q(0., 0., 0.);
      me_K.set_table(*ts, delta);

      double k1 = 2. * M_PI / (double)L;
      std::vector<cx_double> coef_Sz;
      for(int z=0; z < 2; ++z){    
	double kz = M_PI * z;
	for(int x=-L/2; x < L/2; ++x){
	  double kx = k1 * x;
	  for(int y=-L/2; y < L/2; ++y){
	    double ky = k1 * y;

	    /* Checking if the wavevector is inside the BZ. */
	    double factor = BZ_factor_square_half_filling(kx, ky);
	    if ( std::abs(factor) < 1e-12 ) { continue; }

	    /* Eigenenergy */
	    cx_double ek1 = ts->ek1(kx, ky, kz);
	    cx_double tz = ts->tz; 
	    cx_double ek23 = ts->ek23(kx, ky, kz);
	    cx_double ekz = ts->ekz(kx, ky, kz);	      
	    double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, tz, kz, delta);
	    double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, tz, kz, delta);		  
  
	    cx_double denom1 = ek_plus - ek_minus - omega_i_shifted;
	    cx_double denom2 = ek_plus - ek_minus + omega_f_shifted;
  
	    /* Getting the matrix elements */
	    cx_double R[2];
	    if ( me_K.is_table_set() ) {
	      me_K.get_elem(kx, ky, kz, R);
	    } else {
	      me_K.calc_mat_elems(*ts, delta, kx, ky, kz, R);
	    }

	    /* Coefficient of the Sz-like term */
	    cx_double rk = R[0] / denom1 + R[1] / denom2;
	    double z_k = zk(ek1, tz, kz, delta);		  
	    cx_double wk = 0.5 * sqrt(1. - z_k * z_k);		  
	    cx_double coef = factor * rk * wk;
	      
	    coef_Sz.push_back(coef);
	  }  /* end for y */
	}  /* end for x */
      }  /* end for z */

      for(int x0=-L/2; x0 <= L/2; ++x0){
	for(int y0=-L/2; y0 <= L/2; ++y0){
	  for(int z0=0; z0 <= 1; ++z0){
	    if ( (x0 + y0 + z0) & 1 ) { continue; }
	      
	    /* Fourier transform */
	    cx_double sum = 0;	    
	    std::size_t index_k = 0;	    
	    for(int z=0; z < 2; ++z){    
	      double kz = M_PI * z;
	      for(int x=-L/2; x < L/2; ++x){
		double kx = k1 * x;
		for(int y=-L/2; y < L/2; ++y){
		  double ky = k1 * y;

		  /* Checking if the wavevector is inside the BZ. */
		  double factor = BZ_factor_square_half_filling(kx, ky);
		  if ( std::abs(factor) < 1e-12 ) { continue; }

		  double inner_prod_k = kx * x0 + ky * y0 + kz * z0;
		  cx_double phase = exp(1i*inner_prod_k);		  
		  sum += coef_Sz[index_k] * phase;		  
		  ++index_k;
		}  /* end for y */
	      }  /* end for x */
	    }  /* end for z */

	    /* Output */
	    out_coef_Raman << x0 << "  " << y0 << "   " << z0 << "   " << std::real(sum) << "   " << std::imag(sum) << std::endl;	      
	  }  /* end for z0 */
	}  /* end for y0 */	    
      }  /* end for x0 */

      /* Closing the output file. */
      out_coef_Raman.close();	
    }   /* end for nu */
  }   /* end for mu */    
}

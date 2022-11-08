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
#include "calc_Raman.h"
#include "particle_hole1.h"

/* Creating an instance */
// PhiDerIntegrandBilayer pdib;
extern WaveFuncIntegrandBilayer wfib;

/* Considered bonds: +x, +y, +x+y, +x-y, +z */
struct BondDelta {
  int x; int y; int z;
  BondDelta(int dir){
    if ( dir == 0 ) { x = 1; y = 0; z = 0; }
    else if ( dir == 1 ) { x = 0; y = 1; z = 0; }
    else if ( dir == 2 ) { x = 0; y = 0; z = 1; }
    else {
      std::cerr << "Not supported dir value." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  BondDelta(int _x, int _y, int _z):x(_x),y(_y),z(_z){}
};

int inner_prod(BondDelta const& d1, BondDelta const& d2){
  return d1.x * d2.x + d1.y * d2.y + d1.z * d2.z;
}

cx_double bond_to_hopping_bilayer(hoppings_bilayer2 const& ts, BondDelta const& b, int g1, int g2, int sigma){
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
    // cx_double rot_phase = exp(cx_double(0,stag_rot_angle));
    // if (g1 == 1 && g2 == 0) {
    // 	rot_phase = std::conj(rot_phase);	
    // }
    // if ( sigma == down_spin ) {
    // 	rot_phase = std::conj(rot_phase);
    // }
      
    if (b.x + b.y + b.z & 1) {
      if ( std::abs(b.x + b.y) == 1 ) {
	return cx_double(ts.t);
	// return rot_phase * cx_double(ts.t);
	  
      } else if ( std::abs(b.z) == 1 ) {
	if (g2 == 0) {
	  if ( sigma == up_spin ) {
	    return ts.tz;
	  } else {
	    return std::conj(ts.tz);
	  }
	  // return rot_phase * ts.tz;	    
	} else {
	  if ( sigma == up_spin ) {
	    return std::conj(ts.tz);
	  } else {
	    return ts.tz;
	  }
	  // return std::conj(rot_phase * ts.tz);
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

cx_double velocity_U1(hoppings_bilayer2 const& ts, cx_mat const& Udg, cx_mat const& U, double kx, double ky, double kz, BondDelta _bond, int g1, int g2, int sign_m, int sign_n, int sigma_m, int sigma_n){
  int m_idx = sign_spin_index(sign_m, sigma_m);
  int n_idx = sign_spin_index(sign_n, sigma_n);		      
  double inner_prod = kx * _bond.x + ky * _bond.y + kz * _bond.z;
  cx_double phase = exp(1i*inner_prod);
		      
  cx_double v = 0;
  /* Assume there is no hopping with spin flip. */
  for(int _spin: {up_spin, down_spin}){
    int g1_idx = sublattice_spin_index(g1, _spin);
    int g2_idx = sublattice_spin_index(g2, _spin);
    
    cx_double t_21 = bond_to_hopping_bilayer(ts, _bond, g2, g1, _spin); // from g1 to g2
    cx_double Umg2 = Udg(m_idx, g2_idx);
    cx_double Ug1n = U(g1_idx, n_idx);      
    v += t_21 * Umg2 * Ug1n * phase;
    
    cx_double t_12 = bond_to_hopping_bilayer(ts, _bond, g1, g2, _spin); // from g2 to g1
    cx_double Umg1 = Udg(m_idx, g1_idx);
    cx_double Ug2n = U(g2_idx, n_idx);  
    v += - t_12 * Umg1 * Ug2n * std::conj(phase);
  }

  // // for check
  // std::cout << "velocity = " << 1i * v << std::endl;
  
  return 1i * v;
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
  /* Assume that mu does not depend on L for integral over continuous k. */
  double ch_gap, mu;
  std::tie(ch_gap, mu) = calc_charge_gap_bilayer( L, *ts, delta );  /* Finite size */  
  // double mu = calc_chemical_potential_bilayer2( base_dir, L, *ts, delta );  /* Finite size */

  /* Output */
  ofstream out_gap;
  out_gap.open(base_dir / "gap.text");
  out_gap << "Charge gap = " << ch_gap << std::endl;
  out_gap << "mu = " << mu << std::endl;
  
  /* Polarization */
  Polarization Pz( L, L, 2, NSUBL );

  auto calc_gaps = [&](double qx, double qy, double qz){
    /* Setting the ordering vector */
    Pz.set_q(qx, qy, qz);
    if ( !continuous_k ) {
      /* The polarizations are calculated in advance. */
      Pz.set_table( *ts, delta );
    }
    
    /* Calculating gaps */
    double omega_T = 0, omega_L = 0, omega_ph = 0;
    bool return_upper = true;
    // bool verbose = true;
    bool verbose = false;    
    std::tie(omega_T, omega_L, omega_ph) = calc_gap_bilayer(L, *ts, mu, U, T, delta, cbp, Pz, continuous_k, return_upper, verbose);
    
    return std::make_tuple(omega_T, omega_L, omega_ph);
  };

  /* Calculating gaps at (pi, pi, pi)*/
  double omega_T = 0, omega_L = 0, omega_ph = 0;
  std::tie(omega_T, omega_L, omega_ph) = calc_gaps(M_PI, M_PI, M_PI);
  
  /* Output */
  out_gap << "omega_T = " << omega_T << std::endl;
  out_gap << "omega_L = " << omega_L << std::endl;
  out_gap << "omega_ph = " << omega_ph << std::endl;
  
  /* Calculating Psi'*/
  // double Psider = 1.0;
  // std::cout << "Setting Psider to one for simplicity." << std::endl;
  double Psider = calc_Psider(L, *ts, mu, omega_L, delta, cbp, continuous_k);

  /* Output */
  out_gap << "Psider = " << Psider << std::endl;

  /* Calculating the minimum bk square */
  double min_bk_sq = calc_min_bk_sq_bilayer(L, *ts);
    
  /* Constants */
  double stag_rot_angle = 12.0 / 180.0 * M_PI;
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
    // {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0}, {0, 0, 1}
    
    {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, -1, 0},  {0, 0, 1}
    
    // // for check
    // , {-1, 0, 0}, {0, -1, 0}
    
    // {-1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, -1, 0},  {0, 0, 1}            
  };

  auto wavefunc_sublattice_correction = [](int g1, int g2, cx_double psi){
    if (g1 == 0) {
      if (g2 == 0) {
	return psi;
      } else {
	return psi;
      }
    } else {
      if (g2 == 0) {
	return std::conj(psi);
      } else {
	return - psi;
      }
    }
  };
  
  /* Output of the mean field energy gaps. */
  ofstream out_MF_gap;
  out_MF_gap.open(base_dir / "mean_field_energy_gap.text");    
  out_MF_gap << "# Index  kx  ky  kz  E_plus  E_minus  E_gap" << std::endl;  

  std::vector<double> particle_hole_gaps;
  for(int z=0; z < 2; ++z){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; ++x){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; ++y){
	double ky = k1 * y;
	  
	/* Checking if the wavevector is inside the BZ. */
	double factor = BZ_factor(kx, ky);
	if ( std::abs(factor) < 1e-12 ) { continue; }

	/* Eigenenergy */
	cx_double ek1 = ts->ek1(kx, ky, kz);
	cx_double ek23 = ts->ek23(kx, ky, kz);
	cx_double ekz = ts->ekz(kx, ky, kz);	      
	double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts->tz, kz, delta);
	double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts->tz, kz, delta);
	particle_hole_gaps.push_back(ek_plus - ek_minus);
	
	out_MF_gap << kx << "  " << ky << "  " << kz << "  " << ek_plus << "  " << ek_minus << "  " << ek_plus - ek_minus << std::endl;
      }
    }
  }

  /* Sorted mean-field particle-hole energy gaps */
  std::sort(particle_hole_gaps.begin(), particle_hole_gaps.end());
  ofstream out_eh_gap(base_dir / "particle_hole_gaps.text");
  double e = 0;
  for(double eg: particle_hole_gaps){
    if ( std::abs(eg - e) > 1e-10 ) {
      out_eh_gap << eg << std::endl;
      e = eg;
    }
  }
  out_eh_gap.close();
  
  /* Output of the dispersion calculated by RPA. */  
  ofstream out_dispersion_transverse, out_dispersion_longitudinal;
  out_dispersion_transverse.open(base_dir / "dispersion_transverse.text");
  out_dispersion_longitudinal.open(base_dir / "dispersion_longitudinal.text");
  out_dispersion_transverse << "# Index  kx  ky  kz  gap" << std::endl;
  out_dispersion_longitudinal << "# Index  kx  ky  kz  gap" << std::endl;
  
  /* Output of the wave function. */
  ofstream out_wavefunc_k;
  out_wavefunc_k.open(base_dir / "wavefunction_k.text");
  out_wavefunc_k << "# Index  kx  ky  kz   Re[psi_AA]   Im[psi_AA]   Re[psi_AB]   Im[psi_AB]" << std::endl;
  
  /* Energies. */
  std::vector<double> E_k_plus;
  std::vector<double> E_k_minus;
  std::vector<double> gap_T;
  std::vector<double> gap_L;

  /* Setting the parameters */
  int spin = up_spin;  // Unused  
  int sublattice = 0;  // Unused
  int diff_r[] = {0,0,0};  // Unused
		
  std::size_t idx_k_sigma = 0;
  for(int z=0; z < 2; ++z){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; ++x){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; ++y){
	double ky = k1 * y;
	  
	/* Checking if the wavevector is inside the BZ. */
	double factor = BZ_factor(kx, ky);
	if ( std::abs(factor) < 1e-12 ) { continue; }

	for(int sigma: {up_spin, down_spin}){	
	  /* Calculating the gaps */
	  double gT = 0, gL = 0, gph = 0;
	  std::tie(gT, gL, gph) = calc_gaps(kx, ky, kz);  // Solving the pole equation in RPA.
	  gap_T.push_back(gT);
	  gap_L.push_back(gL);

	  /* Eigenenergy */
	  cx_double ek1 = ts->ek1(kx, ky, kz);
	  cx_double ek23 = ts->ek23(kx, ky, kz);
	  cx_double ekz = ts->ekz(kx, ky, kz);	      
	  double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts->tz, kz, delta);
	  double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts->tz, kz, delta);
	  E_k_plus.push_back(ek_plus);
	  E_k_minus.push_back(ek_minus);

	  /* For checking the wave functions. */
	  wfib.set_parameters(*ts, pr.largeUlimit, pr.largeU_scaling_prefactor, sigma, omega_L, Psider, delta, diff_r, sublattice, min_bk_sq);	  
	  cx_double xki = xk(wfib.spin(), ek1, ts->tz, kz, wfib.delta());
	  double zki = zk(ek1, ts->tz, kz, wfib.delta());
	  cx_double bki = bk(wfib.spin(), ek1, ts->tz, kz);
	  cx_double wavefunc_AA = wfib.integrand(xki, zki, bki, ek_plus, ek_minus, 1., 0);
	  cx_double wavefunc_AB = wfib.integrand(xki, zki, bki, ek_plus, ek_minus, 1., 1);	
		    
	  /* Output */
	  int npw = 12;
	  out_dispersion_transverse << idx_k_sigma << std::setw(npw) << kx << std::setw(npw) << ky << std::setw(npw) << kz << std::setw(pw) << gT << std::endl;
	  out_dispersion_longitudinal << idx_k_sigma << std::setw(npw) << kx << std::setw(npw) << ky << std::setw(npw) << kz << std::setw(pw) << gL << std::endl;

	  out_wavefunc_k << idx_k_sigma << std::setw(npw) << kx << std::setw(npw) << ky << std::setw(npw) << kz << std::setw(pw) << std::real(wavefunc_AA) << std::setw(pw) << std::imag(wavefunc_AA) << std::setw(pw) << std::real(wavefunc_AB) << std::setw(pw) << std::imag(wavefunc_AB) << std::endl;;
	
	  ++idx_k_sigma;
	}
      }
    }
  }

  /* Results */
  double *spec_Raman_xx = new double[n_omegas];
  double *spec_Raman_xy = new double[n_omegas];
  double *spec_Raman_xz = new double[n_omegas];
  double *spec_Raman_yx = new double[n_omegas];
  double *spec_Raman_yy = new double[n_omegas];
  double *spec_Raman_yz = new double[n_omegas];
  double *spec_Raman_zx = new double[n_omegas];
  double *spec_Raman_zy = new double[n_omegas];
  double *spec_Raman_zz = new double[n_omegas];        
  
  auto spec_Raman = [&](int mu, int nu){
    double **spec;
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
    return *spec;
  };
  
  /* Calculating the spectral weights and the energies. */    
  // std::vector<std::vector<cx_double>> weight_plus;
  std::vector<std::vector<cx_double>> weight_minus1, weight_minus2;
  for(int mu=0; mu < 3; ++mu){
    BondDelta e_mu(mu);
    for(int nu=0; nu < 3; ++nu){
      BondDelta e_nu(nu);
      // std::vector<cx_double> weight_plus_comp;
      std::vector<cx_double> weight_minus1_comp, weight_minus2_comp;
  
      idx_k_sigma = 0;
      for(int z=0; z < 2; ++z){    
	double kz = M_PI * z;
	for(int x=-L/2; x < L/2; ++x){
	  double kx = k1 * x;
	  for(int y=-L/2; y < L/2; ++y){
	    double ky = k1 * y;

	    /* Checking if the wavevector is inside the BZ. */
	    double factor = BZ_factor(kx, ky);
	    if ( std::abs(factor) < 1e-12 ) { continue; }

	    cx_double ek1 = ts->ek1(kx, ky, kz);
	    cx_mat Uk = gs_HF(ek1, ts->tz, kz, delta);
	    cx_mat Uk_dg = arma::trans(Uk);

	    cx_double ek1_pi = ts->ek1(kx + M_PI, ky + M_PI, kz + M_PI);
	    cx_mat Uk_pi = gs_HF(ek1_pi, ts->tz, kz + M_PI, delta);
	    cx_mat Uk_pi_dg = arma::trans(Uk_pi);

	    // // for check
	    // std::cout << "Uk" << std::endl;
	    // std::cout << Uk << std::endl;

	    for(int sigma: {up_spin, down_spin}){  // Electron spin in the conduction band labeling the final state
	      cx_double mat_elem1_minus = 0.0;
	      cx_double mat_elem2_minus = 0.0;
	      
	      /* ki and kf take 0 and (pi,pi,pi) =: 1. */
	      // std::vector<std::pair<int,int>> ks_photon{{1,0}};
	      // std::vector<std::pair<int,int>> ks_photon{{0,0}};	      
	      std::vector<std::pair<int,int>> ks_photon{{0,0}, {0,1}, {1,0}, {1,1}};
	      
	      for(std::pair<int,int> k_photon: ks_photon){
		int ki = k_photon.first;
		int kf = k_photon.second;
		
		cx_mat *U_ki, *Udg_ki;
		if ( ki == 0 ) {
		  U_ki = &Uk;
		  Udg_ki = &Uk_dg;
		} else {
		  U_ki = &Uk_pi;
		  Udg_ki = &Uk_pi_dg;
		}
		
		cx_mat *U_kf, *Udg_kf;
		if ( kf == 0 ) {
		  U_kf = &Uk;
		  Udg_kf = &Uk_dg;
		} else {
		  U_kf = &Uk_pi;
		  Udg_kf = &Uk_pi_dg;
		}
	      
		/* Velocity */
		cx_double v1 = 0.0, v2 = 0.0;
		for(BondDelta bond1: bonds){
		  int inner_prod_bond1 = inner_prod(e_mu, bond1);
		  if ( inner_prod_bond1 == 0 ) { continue; }
			
		  for(BondDelta bond2: bonds){
		    int inner_prod_bond2 = inner_prod(e_nu, bond2);
		    if ( inner_prod_bond2 == 0 ) { continue; }

		    double dir_comp = (double)(inner_prod_bond1 * inner_prod_bond2);
		  
		    std::vector<std::pair<int,int>> sublattice_pairs{{0,0}, {0,1}, {1,0}, {1,1}};
		    // std::vector<std::pair<int,int>> sublattice_pairs{{0,0}, {0,1}, {1,1}};
		    std::vector<std::pair<int,int>> sublattice_pairs1{{0,1}};
		    std::vector<std::pair<int,int>> sublattice_pairs2{{1,0}};
		    // std::vector<std::pair<int,int>> sublattice_pairs1{{1,0},{0,1}};
		    // std::vector<std::pair<int,int>> sublattice_pairs2{{0,1}, {1,0}};		    
		    
		    // std::vector<std::pair<int,int>> sublattice_pairs{{0,1}};
		    
		    // for(std::pair<int,int> gamma_pair12: sublattice_pairs1){
		    for(std::pair<int,int> gamma_pair12: sublattice_pairs){
		      
		      int g1 = gamma_pair12.first;
		      int g2 = gamma_pair12.second;
		      // for(std::pair<int,int> gamma_pair34: sublattice_pairs2){
		      for(std::pair<int,int> gamma_pair34: sublattice_pairs){
			
			int g3 = gamma_pair34.first;
			int g4 = gamma_pair34.second;
			int sigma_bar = sigma == up_spin ? down_spin : up_spin;			    
			cx_double v_b2_34_pp = velocity_U1(*ts, *Udg_kf, Uk, kx, ky, kz, bond2, g3, g4, 1, 1, sigma, sigma);
			cx_double v_b2_34_mm = velocity_U1(*ts, *Udg_kf, Uk, kx, ky, kz, bond2, g3, g4, -1, -1, sigma_bar, sigma_bar);
			cx_double v_b1_12_pm = velocity_U1(*ts, Uk_dg, *U_ki, kx, ky, kz, bond1, g1, g2,  1, -1, sigma, sigma);
			cx_double v_b1_12_pp = velocity_U1(*ts, *Udg_ki, Uk, kx, ky, kz, bond1, g1, g2,  1,  1, sigma, sigma);
			cx_double v_b1_12_mm = velocity_U1(*ts, *Udg_ki, Uk, kx, ky, kz, bond1, g1, g2, -1, -1, sigma_bar, sigma_bar);
			cx_double v_b2_34_pm = velocity_U1(*ts, Uk_dg, *U_kf, kx, ky, kz, bond2, g3, g4,  1, -1, sigma, sigma);
			// v1 += dir_comp * (v_b2_34_pp + v_b2_34_mm) * v_b1_12_pm;
			// v2 += dir_comp * (v_b1_12_pp + v_b1_12_mm) * v_b2_34_pm;

			cx_double v1i = dir_comp * v_b2_34_pp * v_b1_12_pm;
			cx_double v2i = dir_comp * v_b1_12_pp * v_b2_34_pm;

			// cx_double v1i = dir_comp * (v_b2_34_pp + v_b2_34_mm) * v_b1_12_pm;
			// cx_double v2i = dir_comp * (v_b1_12_pp + v_b1_12_mm) * v_b2_34_pm;
			

			v1 += v1i;
			v2 += v2i;
			
			// // for check
			// std::cout << dir_comp << "  " << g1 << " " << g2 << " " << g3 << " " << g4 << " " << v_b2_34_pp << " " << v_b2_34_mm << " " << v_b1_12_pm << " " << v_b1_12_pp << " " << v_b1_12_mm << " " << v_b2_34_pm << "  " << v1i << " " << v2i << "   " << v1 << " " << v2 << std::endl;

		      }  /* end for gamma_pair34 */			      
		    }  /* end for gamma_pair12 */
		  }  /* end for bond2 */
		}  /* end for bond1 */

		cx_double coef1 = 0.0, coef2 = 0.0;
		for(int g5=0; g5 < 2; ++g5){
		  int g5_idx = sublattice_spin_index(g5, sigma);
		  cx_double U_5plus = Uk(g5_idx, sign_spin_index(1, sigma));
		  
		  for(int g6=0; g6 < 2; ++g6){
		    int g6_idx = sublattice_spin_index(g6, sigma);
		    // cx_double U6_plus = Uk(g6_idx, sign_spin_index(1, sigma));
		    cx_double Udg_minus6 = Uk_dg(sign_spin_index(-1, sigma), g6_idx);

		    
		    /* Using the flipped spin for the wave-function calculation depending on the sublattices. */
		    int sigma_psi = sigma;
		    if ( g6 == 1 && g5 == 0 ) {
		      if ( sigma == up_spin ) {
			sigma_psi = down_spin;
		      } else {
			sigma_psi = up_spin;
		      }
		    }

		    /* Setting the parameters. */
		    int diff_r[] = {0,0,0};  // Unused
		    int sublattice = g5 == g6 ? 0 : 1;		    
		    wfib.set_parameters(*ts, pr.largeUlimit, pr.largeU_scaling_prefactor, sigma_psi, omega_L, Psider, delta, diff_r, sublattice, min_bk_sq);
		
		    /* Calculating Psi(g6, sigma, g5, sigma) from Psi(A, sigma_psi, A/B, sigma_psi) */
		    cx_double xki = xk(wfib.spin(), ek1, ts->tz, kz, wfib.delta());
		    double zki = zk(ek1, ts->tz, kz, wfib.delta());
		    cx_double bki = bk(wfib.spin(), ek1, ts->tz, kz);

		    /* Mean field eigenenergy */
		    double ek_plus = E_k_plus[idx_k_sigma];
		    double ek_minus = E_k_minus[idx_k_sigma];
		  
		    /* Psi(A, sigma_psi, A/B, sigma_psi) */
		    cx_double wavefunc = wfib.integrand(xki, zki, bki, ek_plus, ek_minus, 1., sublattice);

		    // // for check
		    // std::cout << "wavefunc: " << ek_plus << "  " << ek_minus << "  " << xki << "  " << zki << "  " << bki << "  " << sublattice << "  " << wavefunc << std::endl;
		  
		    /* Transforming to Psi(g6, sigma, g5, sigma) */
		    wavefunc = wavefunc_sublattice_correction(g6, g5, wavefunc);

		    /* Adding the matrix elements. */
		    coef1 += U_5plus * Udg_minus6 * std::conj(wavefunc);
		    coef2 += U_5plus * Udg_minus6 * std::conj(wavefunc);		      

		    // // for check
		    // if (std::abs(v1) + std::abs(v2) > 1e-12) {
		    //   std::cout << x << " " << y << " " << z << " " << sigma << " " << g5 << " " << g6 << " " << U5_plus << "  " << U6_minus << "  " << wavefunc << "  " << v1 << "  " << v2 << std::endl;
		    // }
		  
		  }  /* end for g6 */	    
		}  /* end for g5 */

		/* Matrix elements */
		mat_elem1_minus += coef1 * v1;
		mat_elem2_minus += coef2 * v2;
		  
	      }  /* end for k_photon */
	      
	      weight_minus1_comp.push_back(mat_elem1_minus);
	      weight_minus2_comp.push_back(mat_elem2_minus);

	      // // for check
	      // if ( (mu == 0 && nu == 1) || (mu == 1 && nu == 0) ) {
	      // 	// if ( mu == 0 && nu == 1 ) {
		
	      // 	double w1_Re = std::real(mat_elem1_minus);
	      // 	double w1_Im = std::imag(mat_elem1_minus);
	      // 	double w2_Re = std::real(mat_elem2_minus);
	      // 	double w2_Im = std::imag(mat_elem2_minus);				

	      // 	std::cout << mu << "  " << nu << "  " << kx << "  " << ky << "  " << kz << "  " << w1_Re << "  " << w1_Im << "  " << w2_Re << "  " << w2_Im << std::endl;
	      // }
		
	      ++idx_k_sigma;
	    } /* end for sigma */	      
	  } /* end for y */
	} /* end for x */
      } /* end for z */
  
	// weight_plus.push_back(weight_plus_comp);
      weight_minus1.push_back(weight_minus1_comp);
      weight_minus2.push_back(weight_minus2_comp);      
    } /* end for nu */
  } /* end for mu */

  /* Calculating the Raman scattering cross sections. */
  for(int o=0; o < n_omegas; ++o){
    /* Energy of the final photon state */
    double omega_f = pr.omega_i - omegas[o]; // (eV)
    cx_double omega_shifted = cx_double(omegas[o], pr.eta);    
    int idx_comp = 0;    
    for(int mu=0; mu < 3; ++mu){
      for(int nu=0; nu < 3; ++nu){	
	idx_k_sigma = 0;
	cx_double sum = 0.0;	
	for(int z=0; z < 2; ++z){    
	  double kz = M_PI * z;
	  for(int x=-L/2; x < L/2; ++x){
	    double kx = k1 * x;
	    for(int y=-L/2; y < L/2; ++y){
	      double ky = k1 * y;

	      /* Checking if the wavevector is inside the BZ. */
	      double factor = BZ_factor(kx, ky);
	      if ( std::abs(factor) < 1e-12 ) { continue; }

	      for(int sigma: {up_spin, down_spin}){
		/* Mean field eigenenergy */
		double ek_plus = E_k_plus[idx_k_sigma];
		double ek_minus = E_k_minus[idx_k_sigma];		  
		double ek_diff = ek_plus - ek_minus;
		
		cx_double weight1_minus = weight_minus1[idx_comp][idx_k_sigma] / (ek_diff - pr.omega_i);	      
		cx_double weight2_minus = weight_minus2[idx_comp][idx_k_sigma] / (ek_diff + omega_f);
		cx_double weight_minus_sum = weight1_minus + weight2_minus;
		
		// Zero temperature
		double Ei = 0;
		double e_gap = gap_L[idx_k_sigma];	      		
		sum += factor * exp(- beta * Ei) * std::norm(g_photon * weight_minus_sum) / (e_gap - omega_shifted);
	      
		++idx_k_sigma;
	      } /* end for sigma */
	    } /* end for y */
	  } /* end for x */
	} /* end for z */
  
	spec_Raman(mu,nu)[o] = 2.0 * std::imag(sum);
  
	++idx_comp;
      } /* end for nu */
    } /* end for mu */    
  }/* end for o */

  /* Output */
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
  
      ofstream out_raman;
      std::string ofilen("Raman_scattering-"+mu_str+nu_str+".text");
      out_raman.open(base_dir / ofilen);      
      out_raman << "# Omega    R" << std::endl;
  
      /* Output */
      double *spec = spec_Raman(mu,nu);
      for(int o=0; o < n_omegas; ++o){
	double omega = omegas[o];
	out_raman << omega << std::setw(pw) << spec[o] << std::endl;
      } /* end for o */
  
      out_raman.close();        
    } /* end for nu */
  } /* end for mu */

  /* Deleting */
  delete[] spec_Raman_xx;
  delete[] spec_Raman_xy;
  delete[] spec_Raman_xz;
  delete[] spec_Raman_yx;
  delete[] spec_Raman_yy;
  delete[] spec_Raman_yz;
  delete[] spec_Raman_zx;
  delete[] spec_Raman_zy;
  delete[] spec_Raman_zz;
  
  /* Closing the output files. */
  out_gap.close();
  out_dispersion_transverse.close();
  out_dispersion_longitudinal.close();
  out_MF_gap.close();
  out_wavefunc_k.close();
}
			 

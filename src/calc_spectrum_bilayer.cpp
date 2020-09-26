/*****************************************************************************
*
* Functions for calculating the spectrum
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <complex>
#include <omp.h>
#include "calc_spectrum.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"
#include "calc_chemical_potential.h"

std::tuple<double, double> calc_single_particle_energy(hoppings const& ts, double kx, double ky, double kz, double delta){
  double ek1 = ts.ek1(kx, ky, kz);
  double ek2 = ts.ek2(kx, ky, kz);
  double ek3 = ts.ek3(kx, ky, kz);
  double Em = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double Ep = eigenenergy_HF_plus(ek1, ek2, ek3, delta);
  return std::make_tuple(Em, Ep);
}

std::tuple<double, double> calc_single_particle_energy2(hoppings2 const& ts, double kx, double ky, double kz, double delta){
  cx_double ek1 = ts.ek1(kx, ky, kz);
  cx_double ek23 = ts.ek23(kx, ky, kz);
  cx_double ekz = ts.ekz(kx, ky, kz);
  cx_double tz = ts.tz;
  double Em = eigenenergy_HF(-1, ek1, ek23, ekz, tz, kz, delta);
  double Ep = eigenenergy_HF(1, ek1, ek23, ekz, tz, kz, delta);
  return std::make_tuple(Em, Ep);
}

void calc_single_particle_energy_bilayer(hoppings const& ts, int L, double delta){
  /* Output */
  std::ofstream out_e;
  out_e.open("single_particle_energy.text");
  
  /* Wavenumbers */
  double qx = 0;
  double qy = 0;
  double qz = 0;
  int q_idx = 0;

  int prec = 15;
  double k1 = 2. * M_PI / (double)L;
  
  auto output_energy = [&](){
    auto [Em, Ep] = calc_single_particle_energy( ts, qx, qy, qz, delta );
    out_e << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << Em << std::setw( prec ) << Ep << std::endl;
  };
  
  for(int z=0; z < 2; z++){
    qz = M_PI * z;
    
    for(int x=0; x < L/4; x++){
      qx = M_PI - k1 * x;
      qy = k1 * x;
      output_energy();
      ++q_idx;
    }
    
    for(int x=0; x < L/4; x++){
      qx = 0.5 * M_PI - k1 * x;
      qy = 0.5 * M_PI - k1 * x;
      output_energy();
      ++q_idx;
    }
    
    for(int x=0; x < L/2; x++){
      qx = k1 * x;
      qy = 0;
      output_energy();
      ++q_idx;
    }
    
    for(int y=0; y < L/2; y++){
      qx = M_PI;
      qy = k1 * y;
      output_energy();
      ++q_idx;
    }
    
    for(int x=0; x <= L/4; x++){
      qx = M_PI - k1 * x;
      qy = M_PI - k1 * x;
      output_energy();
      ++q_idx;
    }
  }

  out_e.close();
}

void calc_single_particle_energy_bilayer2(hoppings2 const& ts, int L, double delta){
  /* Output */
  std::ofstream out_e;
  out_e.open("single_particle_energy.text");
  
  /* Wavenumbers */
  double qx = 0;
  double qy = 0;
  double qz = 0;
  int q_idx = 0;

  int prec = 15;
  double k1 = 2. * M_PI / (double)L;
  
  auto output_energy = [&](){
    auto [Em, Ep] = calc_single_particle_energy2( ts, qx, qy, qz, delta );
    out_e << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << Em << std::setw( prec ) << Ep << std::endl;
  };
  
  for(int z=0; z < 2; z++){
    qz = M_PI * z;
    
    for(int x=0; x < L/4; x++){
      qx = M_PI - k1 * x;
      qy = k1 * x;
      output_energy();
      ++q_idx;
    }
    
    for(int x=0; x < L/4; x++){
      qx = 0.5 * M_PI - k1 * x;
      qy = 0.5 * M_PI - k1 * x;
      output_energy();
      ++q_idx;
    }
    
    for(int x=0; x < L/2; x++){
      qx = k1 * x;
      qy = 0;
      output_energy();
      ++q_idx;
    }
    
    for(int y=0; y < L/2; y++){
      qx = M_PI;
      qy = k1 * y;
      output_energy();
      ++q_idx;
    }
    
    for(int x=0; x <= L/4; x++){
      qx = M_PI - k1 * x;
      qy = M_PI - k1 * x;
      output_energy();
      ++q_idx;
    }
  }

  out_e.close();
}

void calc_spectrum_bilayer(double theta, double phi, double t3, double U, int L, double eta){
  std::unique_ptr<hoppings_bilayer> ts;
  ts = hoppings_Sr3Ir2O7::mk_Sr3Ir2O7(theta, phi, t3);  
  
  double k1 = 2. * M_PI / (double)L;
  int prec = 15;
  
  /* Omegas */
  double delta_omega = 0.001;  
  double max_omega = ts->t_max() * 10;
  
  int n_omegas = int(max_omega/delta_omega+0.5);
  std::vector<double> omegas(n_omegas);
  for(int o=1; o <= n_omegas; o++){ omegas[o-1] = delta_omega * o; }

  /* Calculate the chemical potential and the charge gap. */
  double delta = solve_self_consistent_eq_bilayer( L, *ts, U );
  std::cout << "delta = " << delta << std::endl;
  double mu = calc_chemical_potential_bilayer( L, *ts, delta );  /* Updated */

  /* Single particle energy */
  calc_single_particle_energy_bilayer( *ts, L, delta );
  
  /* Output */
  std::ofstream out_xy, out_z;
  out_xy.open("spectrum-xy.text");
  out_z.open("spectrum-z.text");

  /* Wavenumbers */
  double qx = 0;
  double qy = 0;
  double qz = 0;
  int q_idx = 0;
  
  /* 2-layer Square lattice */
  /* X -> Sigma -> Gamma -> X -> M -> Sigma */
  /* Gamma = ( 0, 0 )               */
  /* X = ( pi, 0 )                  */
  /* Sigma = ( pi/2, pi/2 )         */
  /* M = ( pi, pi )                 */

  /* From the response function to the dynamic structure factor */  
  auto output_spectrum = [&](){
    std::cout << "( qidx, qx, qy, qz ) = ( " << q_idx << ", " << qx << ", " << qy << ", " << qz << " )" << std::endl;
    for(int o=0; o < n_omegas; o++){
      cx_double cx_omega(omegas[o], eta);

      /* The (S^+ S^-) response function, or the retarded Green's function */
      double factor_dsf = 2.;
      cx_double chi_xy = calc_intensity_bilayer( L, *ts, mu, U, delta, qx, qy, qz, cx_omega, false );
      
      double spec_xy = factor_dsf * std::imag(chi_xy);
      
      /* Output */
      out_xy << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << omegas[o] << std::setw( prec ) << spec_xy << std::setw( prec ) << U << std::endl;

      /* The (S^z S^z) response function, or the retarded Green's function */
      cx_double chi_z = calc_intensity_bilayer( L, *ts, mu, U, delta, qx, qy, qz, cx_omega, true );
      
      double spec_z = factor_dsf * std::imag(chi_z);
      
      /* Output */
      out_z << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << omegas[o] << std::setw( prec ) << spec_z << std::setw( prec ) << U << std::endl;      
    }
  };

  /* Through symmetric points */
  for(int z=0; z < 2; z++){
    qz = M_PI * z;
  
    for(int x=0; x < L/4; x++){
      qx = M_PI - k1 * x;
      qy = k1 * x;
      output_spectrum();
      ++q_idx;
    }
  
    for(int x=0; x < L/4; x++){
      qx = 0.5 * M_PI - k1 * x;
      qy = 0.5 * M_PI - k1 * x;
      output_spectrum();
      ++q_idx;
    }
  
    for(int x=0; x < L/2; x++){
      qx = k1 * x;
      qy = 0;
      output_spectrum();
      ++q_idx;
    }
  
    for(int y=0; y < L/2; y++){
      qx = M_PI;
      qy = k1 * y;
      output_spectrum();
      ++q_idx;
    }
  
    for(int x=0; x <= L/4; x++){
      qx = M_PI - k1 * x;
      qy = M_PI - k1 * x;
      output_spectrum();
      ++q_idx;
    }
  }

  out_xy.close();
  out_z.close();
}

class WaveVector {
public:
  explicit WaveVector(std::string const& _q_type, int _qi, int _qf, int _L){
    q_type_ = _q_type;
    L_ = _L;
    make_q_table();
    
    qi_ = _qi;
    qf_ = _qf < 0 || _qf >= qx_.size() ? qx_.size() - 1 : _qf;
    q_idx_ = qi_;
  }
  std::string q_type() const { return q_type_; }
  int q_idx() const { return q_idx_; }
  int qf() const { return qf_; }
  void increment_index() { ++q_idx_; }
  bool is_done() { return q_idx() > qf() ? true : false; }
  int L() const { return L_; }
  double qx(int i) const { return qx_[i]; }
  double qy(int i) const { return qy_[i]; }
  double qz(int i) const { return qz_[i]; }
  double qx() const { return qx(q_idx()); }
  double qy() const { return qy(q_idx()); }
  double qz() const { return qz(q_idx()); }
  void push_back_qs(double qx, double qy, double qz);
  void make_q_table();
    
private:
  std::string q_type_;
  int q_idx_;
  int qi_;
  int qf_;
  int L_;
  std::vector<double> qx_;
  std::vector<double> qy_;
  std::vector<double> qz_;
};

void WaveVector::push_back_qs(double qx, double qy, double qz){
  qx_.push_back(qx);
  qy_.push_back(qy);
  qz_.push_back(qz);
}

void WaveVector::make_q_table(){
  /* 2-layer Square lattice */
  /* X -> Sigma -> Gamma -> X -> M -> Sigma */
  /* Gamma = ( 0, 0 )               */
  /* X = ( pi, 0 )                  */
  /* Sigma = ( pi/2, pi/2 )         */
  /* M = ( pi, pi )                 */
  
  double qx, qy, qz;
  qx_.clear();
  qy_.clear();
  qz_.clear();
  
  if ( q_type() == "bilayer" ) {
    /* Delta k */
    double k1 = 2. * M_PI / (double)L();

    /* Through symmetric points */  
    for(int z=0; z < 2; z++){
      qz = M_PI * z;
    
      for(int x=0; x < L()/4; x++){
	qx = M_PI - k1 * x;
	qy = k1 * x;
	push_back_qs(qx, qy, qz);
      }
  
      for(int x=0; x < L()/4; x++){
	qx = 0.5 * M_PI - k1 * x;
	qy = 0.5 * M_PI - k1 * x;
	push_back_qs(qx, qy, qz);
      }
  
      for(int x=0; x < L()/2; x++){
	qx = k1 * x;
	qy = 0;
	push_back_qs(qx, qy, qz);
      }
  
      for(int y=0; y < L()/2; y++){
	qx = M_PI;
	qy = k1 * y;
	push_back_qs(qx, qy, qz);      
      }
  
      for(int x=0; x <= L()/4; x++){
	qx = M_PI - k1 * x;
	qy = M_PI - k1 * x;
	push_back_qs(qx, qy, qz);      
      }
    }
  } else {
    std::cerr << "\"q_type\" == " << q_type() << " is not supported.\n";
    std::exit(EXIT_FAILURE);
  }
}

void calc_spectrum_bilayer2(rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  int Lk = pr.Lk;
  double eta = pr.eta;
  double U = pr.U;
  bool continuous_k = pr.continuous_k;

  /* Wavevectors */
  WaveVector wv("bilayer", pr.qi, pr.qf, pr.Lk);
  
  /* Hopping parameters */
  double t1 = pr.t1;
  double t2 = pr.t2;
  double t3 = pr.t3;
  double tz1 = pr.t4;
  double tz2 = pr.t5;
  double phase1 = pr.phase1;
  double phasez = pr.phase4;
  
  using namespace std::complex_literals;  
  // cx_double t1_cx( pr.t1, pr.t1_bar );
  // cx_double tz_cx( pr.t4, pr.t4_bar );
  cx_double t1_cx = t1*exp(1i*phase1*0.5);
  cx_double tz1_cx = tz1*exp(1i*phasez*0.5);

  /* Hopping class */
  std::unique_ptr<hoppings_bilayer2> ts;
  ts = hoppings_bilayer2::mk_bilayer2(t1_cx, t2, t3, tz1_cx, tz2);
  
  /* Omegas */
  double omega_min = pr.omega_min;
  double omega_max = pr.omega_max;
  double omega_delta = pr.omega_delta;  
  
  int n_omegas = int((omega_max - omega_min)/omega_delta+0.5);
  std::vector<double> omegas(n_omegas);
  for(int o=1; o <= n_omegas; o++){ omegas[o-1] = omega_min + omega_delta * o; }

  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  /* Calculate the chemical potential and the charge gap. */
  double delta = solve_self_consistent_eq_bilayer2( L, *ts, U, cbp, continuous_k );  
  std::cout << "delta = " << delta << std::endl;  
  /* Assume that mu does not depend on L for integral over continuous k. */  
  double mu = calc_chemical_potential_bilayer2( L, *ts, delta );  /* Finite size */

  /* Single particle energy */
  calc_single_particle_energy_bilayer2( *ts, L, delta );
  
  /* Output */
  std::ofstream out_xy, out_z;
  out_xy.open("spectrum-xy.text");
  out_z.open("spectrum-z.text");

  /* Precision */
  int prec = 15;

  /* Polarization */
  Polarization Pz( L, L, 2, NSUBL );
    
  /* Results */
  double *spec_xy = new double[n_omegas];
  double *spec_z = new double[n_omegas];

  while( !wv.is_done() ){
    /* Extracting the wavevector */
    int q_idx = wv.q_idx();
    double qx = wv.qx();
    double qy = wv.qy();
    double qz = wv.qz();
    wv.increment_index();

    /* Output */
    std::cout << "( qidx, qx, qy, qz ) = ( " << q_idx << ", " << qx << ", " << qy << ", " << qz << " )" << std::endl;    
  
    /* Setting wavenumbers */
    Pz.set_q( qx, qy, qz );

    /* A factor from the response function to the dynamical structure factor */
    double factor_dsf = 2.;

    if ( continuous_k ) {
      /* Single thread here, but Cuba functions will run in parallel automatically. */
      for(int o=0; o < n_omegas; o++){
	/* Calculating the response (Green's) functions */
	cx_double chi_xy, chi_z;	
	std::tie(chi_xy, chi_z) = calc_intensity_bilayer2( L, *ts, mu, U, delta, cbp, Pz, omegas[o]+1i*eta, continuous_k);
	spec_xy[o] = factor_dsf * std::imag(chi_xy);
	spec_z[o] = factor_dsf * std::imag(chi_z);	
      }
    } else {
      /* The polarizations are calculated in advance. */
      Pz.set_table( *ts, delta );
    
#ifdef WITH_OpenMP
#pragma omp parallel
    {
      /* Assignment for each thread */
      int n_threads = omp_get_num_threads();
      int nt = n_omegas / n_threads;
      int rem = n_omegas % n_threads;      
      int thread_id = omp_get_thread_num();
      int oi, of;
      if ( thread_id < rem ) { nt += 1; }

      /* Results for each thread */
      double *spec_xy_thread = new double[nt];
      double *spec_z_thread = new double[nt];
      
      for(int oidx=0; oidx < nt; oidx++){
      	int o = thread_id + oidx * n_threads; // stride: n_threads
#else
      for(int o=0; o < n_omegas; o++){
#endif
	/* Calculating the response (Green's) functions */
	cx_double chi_xy, chi_z;	
	std::tie(chi_xy, chi_z) = calc_intensity_bilayer2( L, *ts, mu, U, delta, cbp, Pz, omegas[o]+1i*eta, continuous_k);
	
#ifdef WITH_OpenMP
	spec_xy_thread[oidx] = factor_dsf * std::imag(chi_xy);
	spec_z_thread[oidx] = factor_dsf * std::imag(chi_z);
#else
	spec_xy[o] = factor_dsf * std::imag(chi_xy);
	spec_z[o] = factor_dsf * std::imag(chi_z);
#endif
      }
	
#ifdef WITH_OpenMP
      /* Extracting the results */
      for(int oidx=0; oidx < nt; oidx++){
      	int o = thread_id + oidx * n_threads;	
      	spec_xy[o] = spec_xy_thread[oidx];
      	spec_z[o] = spec_z_thread[oidx];	
      }
	
      delete[] spec_xy_thread;
      delete[] spec_z_thread;
      }
#endif
    }
    
    /* Output */
    for(int o=0; o < n_omegas; o++){
      out_xy << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << omegas[o] << std::setw( prec ) << spec_xy[o] << std::setw( prec ) << U << std::endl;
      out_z << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << omegas[o] << std::setw( prec ) << spec_z[o] << std::setw( prec ) << U << std::endl;
    }
  };

  delete[] spec_xy;
  delete[] spec_z;

  out_xy.close();
  out_z.close();
}

/*****************************************************************************
*
* Calculation of chi0 with respect to free electrons
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "plot_chi0.h"
#include "rpa_util.h"

cx_double calc_chi0_each(double t, double mu, double kT, double eta, double E, double qx, double qy, double qz, double kx, double ky, double kz){
  double e_k = energy_free_electron( t, mu, kx, ky, kz );
  double e_kpq = energy_free_electron( t, mu, kx + qx, ky + qy, kz + qz );
  double n_k = fermi_density( e_k, kT, mu );
  double n_kpq = fermi_density( e_kpq, kT, mu );
  using namespace std::complex_literals;
  return ( n_kpq - n_k ) / ( E - ( e_kpq - e_k ) + 1i * eta );
}

cx_double calc_chi0_each(double t, double mu, double kT, double eta, double E, double qx, double qy, double kx, double ky){
  double e_k = energy_free_electron( t, mu, kx, ky );
  double e_kpq = energy_free_electron( t, mu, kx + qx, ky + qy );
  double n_k = fermi_density( e_k, kT, mu );
  double n_kpq = fermi_density( e_kpq, kT, mu );
  using namespace std::complex_literals;
  return ( n_kpq - n_k ) / ( E - ( e_kpq - e_k ) + 1i * eta );
}

cx_double calc_chi0(int L, double t, double mu, double kT, double eta, double E, double qx, double qy, double qz){
  double k1 = 2. * M_PI / (double)L;
  cx_double chi0 = 0.;
  for(int x2=-L/2; x2 < L/2; x2++){
    double kx = k1 * x2;
    for(int y2=-L/2; y2 < L/2; y2++){
      double ky = k1 * y2;
      for(int z2=-L/2; z2 < L/2; z2++){
	double kz = k1 * z2;
	chi0 += calc_chi0_each( t, mu, kT, eta, E, qx, qy, qz, kx, ky, kz );
      }
    }
  }
  return chi0 / std::pow( L, 3 );
}

cx_double calc_chi0(int L, double t, double mu, double kT, double eta, double E, double qx, double qy){
  double k1 = 2. * M_PI / (double)L;
  cx_double chi0 = 0;
  for(int x2=-L/2; x2 < L/2; x2++){
    double kx = k1 * x2;
    for(int y2=-L/2; y2 < L/2; y2++){
      double ky = k1 * y2;
      chi0 += calc_chi0_each( t, mu, kT, eta, E, qx, qy, kx, ky );
    }
  }
  return chi0 / std::pow( L, 2 );
}

cx_double RPA_calc(cx_double chi0, double U){
  return chi0 / ( 1. - U * chi0 );
}

void out_calc(boost::filesystem::ofstream& out, int& idx, int L, double t, double U, double mu, double kT, double eta, double E_max, double delta_E, double qx, double qy){
  for(double E = delta_E; E <= E_max; E += delta_E){
    cx_double chi0 = calc_chi0( L, t, mu, kT, eta, E, qx, qy );
    double chi = std::imag( RPA_calc( chi0, U ) );
    double dssf = 2. / ( 1. - exp( - E / kT ) ) * chi;
    int prec = 15;
    out << idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << E << std::setw( prec ) << dssf << std::setw( prec ) << chi << std::endl;
  }
  ++idx;
}

void out_calc(boost::filesystem::ofstream& out, int& idx, int L, double t, double U, double mu, double kT, double eta, double E_max, double delta_E, double qx, double qy, double qz){
  for(double E = delta_E; E <= E_max; E += delta_E){
    cx_double chi0 = calc_chi0( L, t, mu, kT, eta, E, qx, qy, qz );
    double chi = std::imag( RPA_calc( chi0, U ) );
    double dssf = 2. / ( 1. - exp( - E / kT ) ) * chi;
    int prec = 15;
    out << idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << qz << std::setw( prec ) << E << std::setw( prec ) << dssf << std::setw( prec ) << chi << std::endl;
  }
  ++idx;
}

void calc_chi(){
  int L = 8;
  double k1 = 2. * M_PI / L;
  
  double t = 1.;
  double mu = 0;
  double U = 2.0 * t;
  double kT = 0.001 * t;
  
  double eta = 0.02;
  double E_max = 13.0;
  double delta_E = 0.01;
  
  boost::filesystem::ofstream out;
  out.open("dssf-rpa.text");
  
  double qx = 0.;
  double qy = 0.;
  int idx = 0;

#ifdef __SQUARE__
  out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy );  
  for(int x=1; x <= L/2; x++){
    qx = k1 * x;
    out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy );
  }    
  for(int y=1; y <= L/2; y++){
    qy = k1 * y;
    out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy );
  }
#else
  double qz = 0.;
  out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy, qz );  
  for(int x=1; x <= L/2; x++){
    qx = k1 * x;
    out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy, qz );
  }    
  for(int y=1; y <= L/2; y++){
    qy = k1 * y;
    out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy, qz );
  }
  for(int z=1; z <= L/2; z++){
    qz = k1 * z;
    out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy, qz );
  }
  for(int x=L/2-1; x > 0; x--){
    qx = k1 * x;
    qy = k1 * x;
    qz = k1 * x;
    out_calc( out, idx, L, t, U, mu, kT, eta, E_max, delta_E, qx, qy, qz );
  } 
#endif
  out.close();
}

void calc_size_dependence(){  
  double t = 1.;
  double mu = 0;
  double kT = 1e-4 * t;
  double E = 0.08;
  
  boost::filesystem::ofstream out;
  out.open("chi0-L.text");
  
  double qx = 0;
  double qy = 0;
  // double qz = 0;
  
  const int Ls_size = 11;
  int Ls[ Ls_size ] = { 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768 };

  double Leta = 30.0;

  double q_shift = 2. * M_PI / 16.;
    // double q_shift = 2. * M_PI / 1024.;
  qx = M_PI;
  qy = M_PI - q_shift;
  // qz = 0;

  int prec = 15;
  
  for(int Li=0; Li < Ls_size; Li++){
    int L = Ls[ Li ];
    double eta = Leta / (double)L;
    cx_double chi0 = calc_chi0( L, t, mu, kT, eta, E, qx, qy );
    double x = 1. / (double)L;
   
    out << x << std::setw( prec ) << std::imag( chi0 ) << std::endl;
  }

  out.close();
}

void plot_chi0(){
  double t = 1.;
  double mu = 0;
  double kT = 1e-4 * t; // kT needs to be set small enough
  
  boost::filesystem::ofstream out;
  out.open("chi0-E.text");
  
  double qx = 0;
  double qy = 0;
  // double qz = 0;
  
  double Leta = 30.0; // L * eta is set a constant
  double eps_chi = 0.02;

  double q_shift = 2. * M_PI / 16;
  qx = M_PI;
  qy = M_PI - q_shift;
  // qz = 0;

  int prec = 15;
  // double c = sqrt( 2. );
  // double E0 = c * q_shift;
  double delta_E = 0.02;
  double E0 = delta_E;
  // double chi1_pre = 0.;
  for(double E = E0; E < 10.0; E += delta_E){
    double chi1 = 0.;
    double chi2 = 1000.;
    int L = 64;
    
    while( std::abs( chi2 - chi1 ) / chi1 > eps_chi ){
      chi2 = chi1;
      L *= 2;
      double eta = Leta / (double)L;
      chi1 = std::imag( calc_chi0( L, t, mu, kT, eta, E, qx, qy ) ); // for a square lattice

      // for check
      std::cout << E << std::setw( 10 ) << L << std::setw( prec ) << chi1 << std::endl;
    };
    
    // if ( chi1_pre > chi1 ) { break; }

    out << E << std::setw( prec ) << chi1 << std::setw( prec ) << 1. / chi1 << std::endl;
    // chi1_pre = chi1;

  }
  out.close();
}

void plot_chi0(int L, double eta){
  double t = 1.;
  double mu = 0;
  double kT = 1e-4 * t; // kT needs to be set small enough

  int prec = 15;
      
  /* Omegas */
  double delta_omega = 0.01;
  double max_omega = 10.0;
  int n_omegas = int(max_omega/delta_omega+1);
  std::vector<double> omegas(n_omegas);
  for(int o=0; o < n_omegas; o++){ omegas[o] = delta_omega * o; }
 
  boost::filesystem::ofstream out;
  out.open("chi0-E.text");

  double qx = 0;
  double qy = 0;
  int q_idx = 0;
        
  /* From the response function to the dynamic structure factor */
  double factor_dsf = 2.;

  auto output_spectrum = [&](){
    std::cout << "( qx, qy ) = ( " << qx << ", " << qy << " )" << std::endl;
    for(int o=0; o < n_omegas; o++){
      double sqo = factor_dsf * std::imag( calc_chi0( L, t, mu, kT, eta, omegas[o], qx, qy ) );
      out << q_idx << std::setw( prec ) << qx << std::setw( prec ) << qy << std::setw( prec ) << omegas[o] << std::setw( prec ) << sqo << std::endl;
    }
  };

  /* Calculating at symmetric wavevectors */
  double k1 = 2. * M_PI / (double)L;
  
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

  out.close();
}

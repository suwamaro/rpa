/*****************************************************************************
*
* Functions for calculating the spectrum
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

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
  std::cout << "Setting delta = 0." << std::endl;  /* Gap parameter */
  double delta = 0;  /* to get the Fermi energy */
  double mu = calc_chemical_potential_bilayer( L, *ts, delta );  
  delta = solve_self_consistent_eq_bilayer( L, *ts, mu, U );
  std::cout << "delta = " << delta << std::endl;
  mu = calc_chemical_potential_bilayer( L, *ts, delta );  /* Updated */

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

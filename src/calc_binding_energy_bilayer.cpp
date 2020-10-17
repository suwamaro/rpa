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
#include "calc_binding_energy.h"
#include "BinarySearch.h"

void calc_binding_energy_bilayer(path& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  int Lk = pr.Lk;
  double eta = pr.eta;
  bool continuous_k = pr.continuous_k;

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
  
  /* Parameters for Cuba */
  CubaParam cbp(pr);

  /* Precision */
  int prec = 12;
    
  /* Output */
  ofstream out_be;
  out_be.open(base_dir / "binding-energy.text");
  out_be << "# U    ch_gap    omega_ph    omega_T    (omega_ph - omega_T)    omega_L    (omega_ph - omega_L)" << std::endl;
  out_be << std::setprecision(prec);
  
  double Ui = 2.0;
  double Uf = 20.0;
  double U_delta = 2.0;
  
  for(double U = Ui; U <= Uf; U += U_delta){
    /* Calculate the chemical potential and the charge gap. */
    double delta = solve_self_consistent_eq_bilayer2( L, *ts, U, cbp, continuous_k );  
    std::cout << "delta = " << delta << std::endl;  
    /* Assume that mu does not depend on L for integral over continuous k. */
    double ch_gap, mu;
    std::tie(ch_gap, mu) = calc_charge_gap_bilayer( L, *ts, delta );  /* Finite size */

    /* Polarization */
    Polarization Pz( L, L, 2, NSUBL );
  
    /* Setting the ordering vector */
    Pz.set_q( M_PI, M_PI, M_PI );
    if ( !continuous_k ) {
      /* The polarizations are calculated in advance. */
      Pz.set_table( *ts, delta );
    }

    /* Calculating gaps */    
    double omega_T, omega_L, omega_ph;
    std::tie(omega_T, omega_L, omega_ph) = calc_gap_bilayer(L, *ts, mu, U, delta, cbp, Pz, continuous_k);

    /* Output */
    int space = prec + 8;
    out_be << U << std::setw(space) << ch_gap << std::setw(space) << omega_ph << std::setw(space) << omega_T << std::setw(space) << omega_ph - omega_T << std::setw(space) << omega_L << std::setw(space) << omega_ph - omega_L << std::endl;
  }

  out_be.close();
}

/*****************************************************************************
*
* Functions for obtaining the phase boundary.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_phase_boundary.h"
#include "find_critical_point.h"
#include "find_critical_U.h"

void calc_phase_boundary_U_bilayer(path& base_dir, rpa::parameters& pr){
  std::cout << "Obtaining the phase boundary as a function of U..." << std::endl;
    
  /* Output */
  ofstream out_pb;
  out_pb.open( base_dir / "phase-boundary-U.out");
  out_pb << "# U     tz_c" << std::endl;
  
  /* Precision */
  int prec = 12;
  
  /* U values */
  int n_U = int((pr.U_max - pr.U_min) / pr.U_delta + 1e-12) + 1;
  std::vector<double> Us(n_U);
  for(int Ui=0; Ui < n_U; Ui++){ Us[Ui] = pr.U_min + pr.U_delta * Ui; }

  /* For each U */
  for(int Ui=0; Ui < n_U; Ui++){
    double U = Us[Ui];
    if ( std::abs(U) < 1e-12 ) {
      std::cout << "Skipping the case of U = 0." << std::endl;
      continue;
    }
    pr.U = U;
    if ( pr.fix_J ) {
      pr.t1 = 0.5*sqrt(pr.J*U);
    }
    double tz_c = find_critical_point_bilayer(pr);
    out_pb << std::setprecision(prec) << U << std::setw(prec+8) << tz_c << std::endl;
  }
  
  out_pb.close();
}


void calc_phase_boundary_t4_bilayer(path& base_dir, rpa::parameters& pr){
  std::cout << "Obtaining the phase boundary as a function of t4..." << std::endl;
    
  /* Output */
  ofstream out_pb;
  out_pb.open( base_dir / "phase-boundary-t4.out");
  out_pb << "# tz     U_c" << std::endl;
  
  /* Precision */
  int prec = 12;
  int sw = prec + 10;
  
  /* t4 values */
  int n_t4 = int((pr.t4_max - pr.t4_min) / pr.t4_delta + 1e-12) + 1;
  std::vector<double> t4s(n_t4);
  for(int t4i=0; t4i < n_t4; t4i++){ t4s[t4i] = pr.t4_min + pr.t4_delta * t4i; }

  /* For each t4 */
  for(int t4i=0; t4i < n_t4; t4i++){
    double t4 = t4s[t4i];
    if ( std::abs(t4) < 1e-12 ) {
      std::cout << "Skipping the case of t4 = 0." << std::endl;
      continue;
    }

    /* Setting t4 */
    pr.t4 = t4;
    double shift = std::pow(10, 6);    
    double t4r = std::round(t4*shift)/shift;
    
    if (pr.check_mean_field_function) {
      std::string ofn = "mean_field_eq-t4_"+std::to_string(t4r)+".text";
      std::string out_mff(base_dir/ofn);
      check_mean_field_eq_bilayer(out_mff, pr);
    }    

    /* Output */
    std::string ofn_ = "free_energy-t4_"+std::to_string(t4r)+".text";
    std::string ofn(base_dir/ofn_);    
    ofstream out_e(ofn);
    out_e << "# tz     U   F_disorder   F_order   delta_order   mu_order   ch_gap   diff" << std::endl;

    double Uc = find_critical_U_bilayer(pr);    
    std::cout << "Uc = " << Uc << std::endl;
    std::cout << "Checking the energies." << std::endl;
    double U_delta = 0.00001;
    double Ui = Uc + U_delta;
    double U_max = Uc * 1.0002;
    double delta = 0., mu = 0., ch_gap = 0., diff = 0.;
    bool set_init_val = false;
    for(double U = 0.22; U > Uc; U -= 0.005){
    // for(double U = Ui; U <= U_max; U += U_delta){
      double F1 = 0., F2 = 0.;
      std::tie(F1, F2, delta, mu, ch_gap, diff) = calc_total_energies(pr, U, delta, mu, set_init_val);
      set_init_val = true;
      out_e << std::setw(sw) << t4 << std::setw(sw) << U << std::setw(sw) << F1 << std::setw(sw) << F2 << std::setw(sw) << delta << std::setw(sw) << mu << std::setw(sw) << ch_gap << std::setw(sw) << diff << std::endl;
    }
    out_e.close();
    
    out_pb << std::setprecision(prec) << t4 << std::setw(sw) << Uc << std::endl;
  }

  out_pb.close();  
}

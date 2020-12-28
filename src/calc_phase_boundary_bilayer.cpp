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

void calc_phase_boundary_bilayer(path& base_dir, rpa::parameters& pr){
  std::cout << "Obtaining the phase boundary..." << std::endl;
    
  /* Output */
  ofstream out_pb;
  out_pb.open( base_dir / "phase-boundary.text");
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
    if ( U == 0 ) {
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
  
  // /* t4 values */
  // int n_t4 = int((pr.t4_max - pr.t4_min) / pr.t4_delta + 1e-12) + 1;
  // std::vector<double> t4s(n_t4);
  // for(int t4i=0; t4i < n_t4; t4i++){ t4s[t4i] = pr.t4_min + pr.t4_delta * t4i; }

  // /* For each t4 */
  // for(int t4i=0; t4i < n_t4; t4i++){
  //   double t4 = t4s[t4i];
  //   pr.t4 = t4;
  //   double Uc = find_critical_U_bilayer(pr);
  //   out_pb << std::setprecision(prec) << t4 << std::setw(prec+8) << Uc << std::endl;
  // }

  out_pb.close();  
}

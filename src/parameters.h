/*****************************************************************************
*
* Parameters for the RPA.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _RPA_PARAMETERS_
#define _RPA_PARAMETERS_

#include <string>
#include "rpa.h"

namespace rpa {
  class parameters {
  public:
    parameters();
    explicit parameters(std::string const& ifn);
    double calc_T(double Tc) const;

    /* Choose calculations. */
    bool find_critical_U_bilayer;
    bool find_critical_point_bilayer;
    bool find_critical_T_bilayer;
    bool solve_self_consistent_eqs_bilayer_T;
    bool calc_spectrum_bilayer;
    bool calc_wave_func_bilayer;
    bool calc_binding_energy_bilayer;
    bool calc_phase_boundary_U_bilayer;
    bool calc_phase_boundary_t4_bilayer;
    bool check_mean_field_function;
    bool calc_current_bilayer;    
    bool calc_Raman_bilayer;
    bool calc_Raman_bilayer_coefficient;
    bool calc_mean_field_eigenenergy;    
    bool calc_two_site_problem;
    bool check_details;  // For debug
    bool debug_mode;
    
    int L;  // System size
    int Lx, Ly, Lz;
    int Lk;  // Delta q = 2pi / Lk for plot
    std::string wave_vector_type;
    bool relative_temperature;
    double T;  // Temperature (Kelvin)
    double T_over_Tc;  // T is set to be T_over_Tc * Tc if relative_temperature is true.
    bool T_equal_to_0;
    bool Neel_phase;  // Requiring U > Uc
    double filling;  // Electron filling
    
    /* If continous_k == true, L should not matter to the result. */
    bool continuous_k;  // true: Integral over continous k; L = \infty.    

    /* Parameter for solving the self-consistent equations */
    std::size_t max_iter;
    double epsfunc;
    double mod_prefactor;
    bool use_NewtonRaphson;    
    bool use_NelderMead;
    bool use_1d_solver;
    
    double eta;  // Broadening factor
    double U;  // Onsite Coulomb interaction

    /* Hopping amplitudes */
    double t;
    double t1, t1_bar;  // t*_bar is the imaginary part.
    double t2, t2_bar;
    double t3, t3_bar;
    double t4, t4_bar;
    double t5, t5_bar;
    double t6, t6_bar;

    /* Hopping phases */
    double phase;
    double phase1;
    double phase2;
    double phase3;
    double phase4;
    double phase5;
    double phase6;

    /* Temperatures */
    double T_min;
    double T_max;
    double T_delta;
    
    /* Energy scale of the spectrum */
    double omega_min;
    double omega_max;
    double omega_delta;

    /* Wavevector index: from qi to qf */
    int qi;
    int qf;
    
    /* Parameters for Cuba */
    double epsrel;
    double epsabs;
    int flags;
    int maxeval;
    int key;

    /* Parameters for the calculation of the wavefunction */
    bool largeUlimit;
    double largeU_scaling_prefactor;

    /* Parameters for obtaining the U-tz phase diagram */
    bool fix_J;
    double J;
    double t4_min, t4_max, t4_delta;
    double U_min, U_max, U_delta;
    double init_value;
    bool find_metal_insulator_transition;

    /* Parameters for finding a 1st-order transition point */
    bool find_first_order_transition;    
    bool find_U1st_anneal;
    double find_U1st_U_max;
    double find_U1st_U_min;    
    double find_U1st_U_delta;
    double find_U1st_delta_max;
    double find_U1st_delta_min;    
    double find_U1st_delta_delta;    
    
    /* Parameters for the Raman scattering */
    double omega_i;   // Energy of the initial photon state
    double eta_res;   // Broadening factor for resonant contributions
    double factor_resonant;   // A factor for the resonant contributions
    // int n_ex;   // Maximum number of excitons in the initial state.
    double Omega;   // Energy in the spectrum for the coefficient calculation.
  };

  std::tuple<path, rpa::parameters> extract_parameters(const char* dirn);
}

#endif // _RPA_PARAMETERS_

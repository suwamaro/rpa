/*****************************************************************************
*
* Parameters for the RPA.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <iostream>
#include "parameters.h"
#include "cpptoml.h"

namespace rpa {
  parameters::parameters(){};
  parameters::parameters(std::string const& ifn){
    auto config = cpptoml::parse_file(ifn);
    
    find_critical_U_bilayer = config->get_as<bool>("find_critical_U_bilayer").value_or(false);
    find_critical_point_bilayer = config->get_as<bool>("find_critical_point_bilayer").value_or(false);
    find_critical_T_bilayer = config->get_as<bool>("find_critical_T_bilayer").value_or(false);
    solve_self_consistent_eqs_bilayer_T = config->get_as<bool>("solve_self_consistent_eqs_bilayer_T").value_or(false);
    calc_spectrum_bilayer = config->get_as<bool>("calc_spectrum_bilayer").value_or(false);
    calc_wave_func_bilayer = config->get_as<bool>("calc_wave_func_bilayer").value_or(false);
    calc_binding_energy_bilayer = config->get_as<bool>("calc_binding_energy_bilayer").value_or(false);
    calc_phase_boundary_U_bilayer = config->get_as<bool>("calc_phase_boundary_U_bilayer").value_or(false);
    calc_phase_boundary_t4_bilayer = config->get_as<bool>("calc_phase_boundary_t4_bilayer").value_or(false);
    check_mean_field_function = config->get_as<bool>("check_mean_field_function").value_or(false);    
    calc_current_bilayer = config->get_as<bool>("calc_current_bilayer").value_or(false);    
    calc_Raman_bilayer = config->get_as<bool>("calc_Raman_bilayer").value_or(false);
    calc_Raman_bilayer_coefficient = config->get_as<bool>("calc_Raman_bilayer_coefficient").value_or(false);
    calc_mean_field_eigenenergy = config->get_as<bool>("calc_mean_field_eigenenergy").value_or(false);    
    calc_two_site_problem = config->get_as<bool>("calc_two_site_problem").value_or(false);
    check_details = config->get_as<bool>("check_details").value_or(false);                

    L = config->get_as<int>("L").value_or(16);
    Lx = config->get_as<int>("Lx").value_or(4);
    Ly = config->get_as<int>("Ly").value_or(4);
    Lz = config->get_as<int>("Lz").value_or(4);    
    Lk = config->get_as<int>("Lk").value_or(L);    
    wave_vector_type = config->get_as<std::string>("wave_vector_type").value_or("high_symmetry1");
    relative_temperature = config->get_as<bool>("relative_temperature").value_or(false);
    if ( relative_temperature ) {
      T_over_Tc = config->get_as<double>("T_over_Tc").value_or(0.8);
      if ( T_over_Tc < 0 ) {
	std::cerr << "Input parameter \"T_over_Tc\" has to be non-negative." << std::endl;
	std::exit(EXIT_FAILURE);
      } else if ( T_over_Tc < 1e-15 ) {
	T_equal_to_0 = true;
      } else {
	T_equal_to_0 = false;
      }
      
      if ( config->contains("T") ) {
	std::cerr << "Input parameter \"T\" will be ignored, and \"T_over_Tc\" will be used instead." << std::endl;
      }
    } else {
      T = config->get_as<double>("T").value_or(0.0);
      if ( T < 0 ) {
	std::cerr << "Input parameter \"T\" has to be non-negative." << std::endl;
	std::exit(EXIT_FAILURE);
      } else if ( T < 1e-15 ) {
	T_equal_to_0 = true;
      } else {
	T_equal_to_0 = false;
      }
    }

    Neel_phase = config->get_as<bool>("Neel_phase").value_or(false);    
    filling = config->get_as<double>("filling").value_or(0.5);    
    continuous_k = config->get_as<bool>("continuous_k").value_or(false);
    max_iter = config->get_as<std::size_t>("max_iter").value_or(200);    
    epsfunc = config->get_as<double>("epsfunc").value_or(1e-10);
    mod_prefactor = config->get_as<double>("mod_prefactor").value_or(1.0);
    use_NelderMead = config->get_as<bool>("use_NelderMead").value_or(false);
    use_1d_solver = config->get_as<bool>("use_1d_solver").value_or(true);        
    eta = config->get_as<double>("eta").value_or(0.001);
    U = config->get_as<double>("U").value_or(1.0);

    /* Hopping amplitudes */
    t = config->get_as<double>("t").value_or(0);    
    t1 = config->get_as<double>("t1").value_or(0);
    t1_bar = config->get_as<double>("t1_bar").value_or(0);
    t2 = config->get_as<double>("t2").value_or(0);
    t2_bar = config->get_as<double>("t2_bar").value_or(0);
    t3 = config->get_as<double>("t3").value_or(0);
    t3_bar = config->get_as<double>("t3_bar").value_or(0);
    t4 = config->get_as<double>("t4").value_or(0);
    t4_bar = config->get_as<double>("t4_bar").value_or(0);
    t5 = config->get_as<double>("t5").value_or(0);
    t5_bar = config->get_as<double>("t5_bar").value_or(0);
    t6 = config->get_as<double>("t6").value_or(0);
    t6_bar = config->get_as<double>("t6_bar").value_or(0);

    /* Hopping phases */
    phase = config->get_as<double>("phase").value_or(0);    
    phase1 = config->get_as<double>("phase1").value_or(0);
    phase2 = config->get_as<double>("phase2").value_or(0);
    phase3 = config->get_as<double>("phase3").value_or(0);
    phase4 = config->get_as<double>("phase4").value_or(0);
    phase5 = config->get_as<double>("phase5").value_or(0);
    phase6 = config->get_as<double>("phase6").value_or(0);

    /* Temperatures */
    T_min = config->get_as<double>("T_min").value_or(0);
    T_max = config->get_as<double>("T_max").value_or(2000.);
    T_delta = config->get_as<double>("T_delta").value_or(100.);
    
    /* Energy scale of the spectrum */
    omega_min = config->get_as<double>("omega_min").value_or(0);
    omega_max = config->get_as<double>("omega_max").value_or(1.);
    omega_delta = config->get_as<double>("omega_delta").value_or(0.01);

    /* Wavevector index: from qi to qf */
    qi = config->get_as<int>("qi").value_or(0);
    qf = config->get_as<int>("qf").value_or(-1);
    
    /* Parameters for Cuba */
    epsrel = config->get_as<double>("epsrel").value_or(1e-3);
    epsabs = config->get_as<double>("epsabs").value_or(1e-10);
    flags = config->get_as<int64_t>("flags").value_or(0);
    maxeval = config->get_as<int64_t>("maxeval").value_or(1<<25);
    key = config->get_as<int64_t>("key").value_or(0);

    /* Parameters for the calculation of the wavefunction */
    largeUlimit = config->get_as<bool>("largeUlimit").value_or(false);
    largeU_scaling_prefactor = config->get_as<double>("largeU_scaling_prefactor").value_or(0.);

    /* Parameters for obtaining the U-tz phase diagram */
    fix_J = config->get_as<bool>("fix_J").value_or(false);
    J = config->get_as<double>("J").value_or(1.0);
    t4_min = config->get_as<double>("t4_min").value_or(0);
    t4_max = config->get_as<double>("t4_max").value_or(1.0);
    t4_delta = config->get_as<double>("t4_delta").value_or(0.1);
    U_min = config->get_as<double>("U_min").value_or(0);
    U_max = config->get_as<double>("U_max").value_or(1.0);
    U_delta = config->get_as<double>("U_delta").value_or(0.1);
    init_value = config->get_as<double>("init_value").value_or(0.1);
    find_metal_insulator_transition = config->get_as<bool>("find_metal_insulator_transition").value_or(false);
    
    /* Parameters for finding a 1st-order transition point */
    find_first_order_transition = config->get_as<bool>("find_first_order_transition").value_or(false);    
    find_U1st_anneal = config->get_as<bool>("find_U1st_anneal").value_or(true);
    find_U1st_U_max = config->get_as<double>("find_U1st_U_max").value_or(0.3);
    find_U1st_U_delta = config->get_as<double>("find_U1st_U_delta").value_or(0.01);        
    
    /* Parameters for the Raman scattering */
    if (calc_Raman_bilayer) {
      double eV_532_inv_nm = planck_h * c_light / (532.0 * 1e-9);  // (eV)    
      omega_i = config->get_as<double>("omega_i").value_or(eV_532_inv_nm);  // (eV)
      eta_res = config->get_as<double>("eta_res").value_or(0.5*eta);
      
      std::cout << "Incident photon energy is " << omega_i << " eV." << std::endl;    
      std::cout << "1 eV corresponds to " << 1e-2 / (planck_h * c_light) << " cm^{-1}." << std::endl;
      
      factor_resonant = config->get_as<double>("factor_resonant").value_or(1.0);  // (eV)    
      // n_ex = config->get_as<int>("n_ex").value_or(1);
    }

    if (calc_Raman_bilayer_coefficient) {    
      Omega = config->get_as<double>("Omega").value_or(0.1);  // (nm)
    }
  }
  
  double parameters::calc_T(double Tc) const {
    if ( relative_temperature ) {
      return T_over_Tc * Tc;
    } else {
      return T;
    }
  }

  std::tuple<path, rpa::parameters> extract_parameters(const char* dirn){
    /* Simulation directory */
    path base_dir(dirn);
  
    /* Input parameters */  
    auto input_file = base_dir / "config.toml";
    if ( !exists(input_file) ) {
      std::cerr << "Required input file " << input_file << " does not exist.\n";
      std::exit(EXIT_FAILURE);
    }  
    rpa::parameters p(input_file.string());
    return std::make_tuple(base_dir, p);
  }
  
}

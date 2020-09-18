/*****************************************************************************
*
* Parameters for the RPA.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "parameters.h"
#include "cpptoml.h"

namespace rpa {
  parameters::parameters(const char* ifn){
    auto config = cpptoml::parse_file(ifn);

    L = config->get_as<int64_t>("L").value_or(16);
    Lk = config->get_as<int64_t>("Lk").value_or(L);
    q_idx = config->get_as<int64_t>("q_idx").value_or(0);
    continuous_k = config->get_as<bool>("continuous_k").value_or(false);
    eta = config->get_as<double>("eta").value_or(0.001);
    U = config->get_as<double>("U").value_or(1.0);

    /* Hopping amplitudes */    
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
    phase1 = config->get_as<double>("phase1").value_or(0);
    phase2 = config->get_as<double>("phase2").value_or(0);
    phase3 = config->get_as<double>("phase3").value_or(0);
    phase4 = config->get_as<double>("phase4").value_or(0);
    phase5 = config->get_as<double>("phase5").value_or(0);
    phase6 = config->get_as<double>("phase6").value_or(0);

    /* Energy scale of the spectrum */
    omega_min = config->get_as<double>("omega_min").value_or(0);
    omega_max = config->get_as<double>("omega_max").value_or(1.);
    omega_delta = config->get_as<double>("omega_delta").value_or(0.01);

    /* Parameters for Cuba */
    epsrel = config->get_as<double>("epsrel").value_or(1e-3);
    epsabs = config->get_as<double>("epsabs").value_or(1e-10);
    flags = config->get_as<int64_t>("flags").value_or(0);
    maxeval = config->get_as<int64_t>("maxeval").value_or(1<<25);
    key = config->get_as<int64_t>("key").value_or(0);        
  }
}

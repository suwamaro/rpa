/*****************************************************************************
*
* Parameters for the RPA.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef __parameters__
#define __parameters__

namespace rpa {
  class parameters {
  public:
    explicit parameters(const char* ifn);
    int L;  // System size
    int Lk;  // Delta q = 2pi / Lk for plot
    
    /* If continous_k == true, L should not matter to the result. */
    bool continuous_k;  // true: Integral over continous k; L = \infty.    
    
    double eta;  // Broadnening factor
    double U;  // Onsite Coulomb interaction

    /* Hopping amplitudes */
    double t1, t1_bar;  // t*_bar is the imaginary part.
    double t2, t2_bar;
    double t3, t3_bar;
    double t4, t4_bar;
    double t5, t5_bar;
    double t6, t6_bar;

    /* Hopping phases */    
    double phase1;
    double phase2;
    double phase3;
    double phase4;
    double phase5;
    double phase6;

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
  };
}

#endif // __parameters__

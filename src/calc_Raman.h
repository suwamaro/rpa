/*****************************************************************************
*
* Functions for calculating the Raman scattering cross section.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _CALC_RAMAN_H
#define _CALC_RAMAN_H

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"

// /* Integrand */
// class PhiDerIntegrand {
// public:
//   explicit PhiDerIntegrand(hoppings2 *ts):ts_(ts){}
//   virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
//   static double integrand(double omega, double delta, double zk, cx_double bk);
//   hoppings2 *ts() const { return ts_; }
//   virtual ~PhiDerIntegrand(){}
  
// private:  
//   hoppings2 *ts_;
// };

// class PhiDerIntegrandBilayer : public PhiDerIntegrand {
// public:
//   PhiDerIntegrandBilayer():PhiDerIntegrand(&hb_){}
//   void set_parameters(hoppings_bilayer2 const& h, double omega, double delta);
//   int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
//   double omega() const { return omega_; }
//   double delta() const { return delta_; }  

// private:
//   hoppings_bilayer2 hb_;
//   double omega_;
//   double delta_;  
// };

// class WaveFuncIntegrand {
// public:
//   explicit WaveFuncIntegrand(hoppings2 *ts):ts_(ts){}
//   virtual int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const = 0;
//   hoppings2 *ts() const { return ts_; }
//   virtual ~WaveFuncIntegrand(){}
  
// private:  
//   hoppings2 *ts_;
// };

// class WaveFuncIntegrandBilayer : public WaveFuncIntegrand {
// public:
//   WaveFuncIntegrandBilayer():WaveFuncIntegrand(&hb_){}
//   void set_parameters(hoppings_bilayer2 const& h, bool largeUlimit, double scaling_prefactor, int spin, double omega, double psider, double delta, int *diff_r, int sublattice, double min_bk_sq);
//   int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const;
//   cx_double integrand(cx_double xk, double zk, cx_double bk, double ek_plus, double ek_minus, cx_double phase) const;
//   bool largeUlimit() const { return largeUlimit_; }
//   double scaling_prefactor() const { return scaling_prefactor_; }
//   int spin() const { return spin_; }
//   double omega() const { return omega_; }
//   double psider() const { return psider_; }
//   double delta() const { return delta_; }
//   int diff_r(int d) const { return diff_r_[d]; }
//   int sublattice() const { return sublattice_; }
//   double min_bk_sq() const { return min_bk_sq_; }

// private:
//   hoppings_bilayer2 hb_;
//   bool largeUlimit_;
//   double scaling_prefactor_;
//   int spin_;
//   double omega_;
//   double psider_;
//   double delta_;
//   int diff_r_[3];  // x, y, and z
//   int sublattice_;
//   double min_bk_sq_;
// };

// /* For bilayer lattices */
// double pole_eq_bilayer(int L, hoppings_bilayer2 const& ts, double omega, double mu, double U, double T, double delta, CubaParam const& cbp, Polarization const& Pz, bool continuous_k, std::string const& mode);
// double solve_pole_eq_bilayer(int L, hoppings_bilayer2 const& ts, double mu, double U, double T, double delta, CubaParam const& cbp, Polarization const& Pz, bool continuous_k, std::string const& mode, double upper);
// std::tuple<double, double, double> calc_gap_bilayer(int L, hoppings_bilayer2 const& ts, double mu, double U, double T, double delta, CubaParam const& cbp, Polarization const& Pz, bool continuous_k);

// void calc_wave_func_bilayer(path& base_dir, rpa::parameters const& pr);
// void check_wave_func_bilayer(path& base_dir, rpa::parameters const& pr);

void calc_Raman_bilayer(path& base_dir, rpa::parameters const& pr);

#endif // _CALC_RAMAN_H

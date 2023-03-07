/*****************************************************************************
*
* RPA calculation for the fermionic Hubbard model
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _RPA_UTIL_H
#define _RPA_UTIL_H

#include "rpa.h"
#include "hoppings.h"

/* Considered bonds: +x, +y, +x+y, +x-y, +z */
struct BondDelta {
  int x; int y; int z;
  BondDelta():x(0),y(0),z(0){}
  explicit BondDelta(int dir){
    if ( dir == 0 ) { x = 1; y = 0; z = 0; }
    else if ( dir == 1 ) { x = 0; y = 1; z = 0; }
    else if ( dir == 2 ) { x = 0; y = 0; z = 1; }
    else {
      std::cerr << "Not supported dir value." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  BondDelta(int _x, int _y, int _z):x(_x),y(_y),z(_z){}
  bool operator ==(BondDelta const& bd) const {
    if ( this->x == bd.x && this->y == bd.y && this->z == bd.z ) { return true; }
    else { return false; }
  }
};

/* Operators */
struct OperatorK {
  OperatorK():Gamma_1plus(NSUBL*NSUBL, NSUBL*NSUBL),
	      Gamma_1minus(NSUBL*NSUBL, NSUBL*NSUBL),	      
	      Gamma_1z(NSUBL*NSUBL, NSUBL*NSUBL),
	      Gamma_2plus(NSUBL*NSUBL, NSUBL*NSUBL),
	      Gamma_2minus(NSUBL*NSUBL, NSUBL*NSUBL),
	      Gamma_2z(NSUBL*NSUBL, NSUBL*NSUBL){
    Gamma_1plus(0,1) = 1.;
    Gamma_1minus(1,0) = 1.;
    Gamma_1z(0,0) = 1.;
    Gamma_1z(1,1) = -1.;
    Gamma_2plus(2,3) = 1.;
    Gamma_2minus(3,2) = 1.;    
    Gamma_2z(2,2) = 1.;
    Gamma_2z(3,3) = -1.;
  }
  sp_mat Gamma_1plus;
  sp_mat Gamma_1minus;
  sp_mat Gamma_1z;
  sp_mat Gamma_2plus;
  sp_mat Gamma_2minus;
  sp_mat Gamma_2z;
};

double energy_free_electron(double t, double mu, double kx, double ky, double kz);
double energy_free_electron(double t, double mu, double kx, double ky);
// double energy_free_electron_bilayer(hoppings const& ts, double mu, double kx, double ky, double kz);
double energy_free_electron_bilayer1(hoppings const& ts, double mu, double kx, double ky, double kz);
double energy_free_electron_bilayer2(hoppings const& ts, double mu, double kx, double ky, double kz);
double eigenenergy_HF_minus(double e_free, double delta);
double eigenenergy_HF_plus(double e_free, double delta);
double eigenenergy_HF_minus(double ek1, double ek2, double ek3, double delta);
double eigenenergy_HF_plus(double ek1, double ek2, double ek3, double delta);
cx_double bk(int spin, cx_double ek1, cx_double tz, double kz);
double bk(cx_double ek1, cx_double tz, double kz);
double zk(int spin, cx_double ek1, cx_double tz, double kz, double delta);
double zk(cx_double ek1, cx_double tz, double kz, double delta);
double zk_over_delta(int spin, cx_double ek1, cx_double tz, double kz, double delta);
double zk_over_delta(cx_double ek1, cx_double tz, double kz, double delta);
cx_double xk(int spin, cx_double ek1, cx_double tz, double kz, double delta);
double eigenenergy_HF(double sign, cx_double ek1, cx_double ek23, cx_double ekz, cx_double tz, double kz, double delta);
double fermi_energy(double x, double kT, double mu);
double fermi_density(double x, double kT, double mu);
double compressibility(double x, double kT, double mu);

#endif // _RPA_UTIL_H

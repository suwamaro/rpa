/*****************************************************************************
*
* Functions for calculating the gap.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_gap.h"
#include "rpa_util.h"
#include "BinarySearch.h"

double calc_band_gap_bilayer(int L, hoppings_bilayer2 const& ts, double delta, double qx, double qy, double qz){
  double k1 = 2. * M_PI / (double)L;
  double E_diff_min = std::numeric_limits<double>::max();
  cx_double tz = ts.tz;
  
  for(int z=0; z < 2; ++z){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; ++x){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; ++y){
	double ky = k1 * y;      

	cx_double ek1 = ts.ek1(kx, ky, kz);
	cx_double ek23 = ts.ek23(kx, ky, kz);
	cx_double ekz = ts.ekz(kx, ky, kz);

	double kx2 = kx + qx;
	double ky2 = ky + qy;
	double kz2 = kz + qz;

	cx_double ek_q1 = ts.ek1(kx2, ky2, kz2);
	cx_double ek_q23 = ts.ek23(kx2, ky2, kz2);
	cx_double ek_qz = ts.ekz(kx2, ky2, kz2);
	
	double Ek_plus = eigenenergy_HF(1., ek_q1, ek_q23, ek_qz, tz, kz2, delta);
	double Ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, tz, kz, delta);

	double E_diff = std::max( 0., Ek_plus - Ek_minus );
	E_diff_min = std::min( E_diff_min, E_diff );
      } /* end for y */
    } /* end for x */
  } /* end for z */
  
  return E_diff_min;
}

double calc_direct_band_gap_bilayer(int L, hoppings_bilayer2 const& ts, double delta, double kx, double ky, double kz){
  cx_double ek1 = ts.ek1(kx, ky, kz);
  cx_double ek23 = ts.ek23(kx, ky, kz);
  cx_double ekz = ts.ekz(kx, ky, kz);
  double Ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts.tz, kz, delta);
  double Ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts.tz, kz, delta);
  double E_diff = std::max( 0., Ek_plus - Ek_minus );  
  return E_diff;
}

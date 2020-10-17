/*****************************************************************************
*
* Functions for calculating the spectrum
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_chemical_potential.h"

double calc_chemical_potential_bilayer(int L, hoppings const& ts, double delta){
  /* Find the chemical potential as the average of the maximum of the lower band and the minimum of the upper. */
  double k1 = 2. * M_PI / (double)L;
  double E_min = std::numeric_limits<double>::max();
  double E_max = - std::numeric_limits<double>::max();

  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;      

	double ek1 = ts.ek1(kx, ky, kz);
	double ek2 = ts.ek2(kx, ky, kz);
	double ek3 = ts.ek3(kx, ky, kz);

	double Ek_up = eigenenergy_HF_plus(ek1, ek2, ek3, delta);
	E_min = std::min( E_min, Ek_up );
	
	double Ek_low = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
	E_max = std::max( E_max, Ek_low );
      } /* end for y */
    } /* end for x */
  } /* end for z */

  double mu = 0.5 * ( E_min + E_max );  
  std::cout << "E_min of the upper band = " << E_min << std::endl;
  std::cout << "E_max of the lower band = " << E_max << std::endl;
  if ( std::abs( E_min - E_max ) < 1e-12 ) {
    std::cout << "The half filling is metallic." << std::endl;
  } else if ( E_min > E_max ) {
    std::cout << "There is a gap: " << E_min - E_max << std::endl;
  } else {
    std::cout << "There is a intersection: " << E_max - E_min << std::endl;
  }
  std::cout << "mu = " << mu << std::endl;  
  return mu;
}

double calc_chemical_potential_bilayer2(path& base_dir, int L, hoppings2 const& ts, double delta){
  /* Find the chemical potential as the average of the maximum of the lower band and the minimum of the upper. */
  double k1 = 2. * M_PI / (double)L;
  double E_min = std::numeric_limits<double>::max();
  double E_max = - std::numeric_limits<double>::max();
  cx_double tz = ts.tz;
  
  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;      

	cx_double ek1 = ts.ek1(kx, ky, kz);
	cx_double ek23 = ts.ek23(kx, ky, kz);
	cx_double ekz = ts.ekz(kx, ky, kz);

	double Ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, tz, kz, delta);
	double Ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, tz, kz, delta);
	
	E_min = std::min( E_min, Ek_plus );
	E_max = std::max( E_max, Ek_minus );
      } /* end for y */
    } /* end for x */
  } /* end for z */

  double mu = 0.5 * ( E_min + E_max );

  /* Output */
  ofstream out_E;
  out_E.open(base_dir / "eigenenergy.text");
  
  out_E << "E_min of the upper band = " << E_min << std::endl;
  out_E << "E_max of the lower band = " << E_max << std::endl;
  if ( std::abs( E_min - E_max ) < 1e-12 ) {
    out_E << "The half filling is metallic." << std::endl;
  } else if ( E_min > E_max ) {
    out_E << "There is a gap: " << E_min - E_max << std::endl;
  } else {
    out_E << "There is a intersection: " << E_max - E_min << std::endl;
  }
  out_E << "mu = " << mu << std::endl;
  
  out_E.close();
  
  return mu;
}

std::tuple<double, double> calc_charge_gap_bilayer(int L, hoppings2 const& ts, double delta){
  /* Find the chemical potential as the average of the maximum of the lower band and the minimum of the upper. */
  double k1 = 2. * M_PI / (double)L;
  double E_min = std::numeric_limits<double>::max();
  double E_max = - std::numeric_limits<double>::max();
  cx_double tz = ts.tz;
  
  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;      

	cx_double ek1 = ts.ek1(kx, ky, kz);
	cx_double ek23 = ts.ek23(kx, ky, kz);
	cx_double ekz = ts.ekz(kx, ky, kz);

	double Ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, tz, kz, delta);
	double Ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, tz, kz, delta);
	
	E_min = std::min( E_min, Ek_plus );
	E_max = std::max( E_max, Ek_minus );
      } /* end for y */
    } /* end for x */
  } /* end for z */

  double ch_gap = std::max( 0., E_min - E_max );  
  double mu = 0.5 * ( E_min + E_max );
  
  return std::make_tuple(ch_gap, mu);
}

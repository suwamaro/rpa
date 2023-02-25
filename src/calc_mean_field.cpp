/*****************************************************************************
*
* Functions for obtaining the mean field solution.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <complex>
#ifdef WITH_OpenMP
#include <omp.h>
#endif
#include "array3.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"
#include "calc_chemical_potential.h"
#include "BinarySearch.h"

void calc_mean_field_eigenenergy(path& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  int L = pr.L;
  double U = pr.U;
  double filling = pr.filling;
  double T = pr.T;
  bool continuous_k = pr.continuous_k;

  /* Hopping parameters for the bilayer system */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);

  /* Parameters for Cuba */
  CubaParam cbp(pr);
  
  /* Calculate the chemical potential and the charge gap. */
  double delta = solve_self_consistent_eq_bilayer2( L, *ts, U, filling, T, cbp, continuous_k );  
  std::cout << "delta = " << delta << std::endl;
  /* Assume that the chemical potential does not depend on L for integral over continuous k. */
  double ch_gap, ch_pot;
  std::tie(ch_gap, ch_pot) = calc_charge_gap_bilayer( L, *ts, delta );  /* Finite size */

  /* Output file */
  ofstream out_energy;
  out_energy.open(base_dir / "mean_field_energy.text");
  out_energy << "# kx   ky   kz   E_plus   E_minus   Gap" << std::endl;
  
  /* Precision */
  int prec = 12;
  int pw = prec + 10;

  std::vector<double> Egap;
  double k1 = 2. * M_PI / (double)L;
  for(int z=0; z < 2; ++z){    
    double kz = M_PI * z;
    for(int y=-L/2; y < L/2; ++y){
      double ky = k1 * y;	
      for(int x=-L/2; x < L/2; ++x){
	double kx = k1 * x;

	/* Checking if the wavevector is inside the BZ. */
	double factor = BZ_factor_square_half_filling(kx, ky);
	if ( std::abs(factor) < 1e-12 ) { continue; }

	/* Eigenenergy */
	cx_double ek1 = ts->ek1(kx, ky, kz);
	cx_double tz = ts->tz; 
	cx_double ek23 = ts->ek23(kx, ky, kz);
	cx_double ekz = ts->ekz(kx, ky, kz);	      
	double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, tz, kz, delta);
	double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, tz, kz, delta);
	double ek_gap = ek_plus - ek_minus;
	Egap.push_back(ek_gap);
	
	/* Output */
	int prec_orig = out_energy.precision();	
	out_energy << std::setw(13) << kx << std::setw(13) << ky << std::setw(13) << kz
		   << std::setprecision(prec)
		   << std::setw(pw) << ek_plus << std::setw(pw) << ek_minus
		   << std::setw(pw) << ek_gap << std::endl
		   << std::setprecision(prec_orig);
      }
    }
  }

  /* Output file */
  ofstream out_gap;
  out_gap.open(base_dir / "mean_field_energy_gap.text");
  out_gap << std::setprecision(prec);

  /* Sorting */
  std::sort(Egap.begin(), Egap.end());

  /* Output */
  for(double Eg: Egap){
    out_gap << Eg << std::endl;
  }

  /* Closing */
  out_energy.close();
  out_gap.close();    
}

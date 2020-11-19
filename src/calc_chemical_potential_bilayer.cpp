/*****************************************************************************
*
* Functions for calculating the spectrum
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "BinarySearch.h"
#include "rpa_util.h"
#include "calc_single_particle_energy.h"
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

/* Creating an instance */
ElecFillingIntegrandBilayer efib;

int ef_integrand_bilayer(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return efib.calc(ndim, xx, ncomp, ff, userdata);
}

double elec_filling_eq_bilayer(int L, hoppings_bilayer2 const& ts, double mu, double T, double delta, CubaParam const& cbp, bool continuous_k){  
  /* Monotonically decreasing as a function of delta */
  double sum = 0;
  
  if ( continuous_k ) {
    /* Changing the parameters */
    efib.set_parameters(ts, mu, T, delta);

    /* For Cuba */
    double epsrel = 1e-10;
    double epsabs = 1e-10;    
    int ncomp = 1;
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, ef_integrand_bilayer, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);    
    // for(int comp = 0; comp < cbp.NCOMP; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);

    sum = integral[0] / 4.;
  } else {
    double eps = 1e-12;
    double k1 = 2. * M_PI / (double)L;
    
    for(int z=0; z < 2; z++){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; x++){
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; y++){
	  double ky = k1 * y;      
	
	  /* Checking if k is inside/outside the BZ. */
	  double mu_free = 0;  /* Assume at half filling */
	  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
	  if ( e_free > mu_free + eps ) continue;

	  /* Prefactor */
	  double factor = 1.;
	
	  /* On the zone boundary */
	  if ( std::abs(e_free - mu_free) < eps ) {	
	    factor = 0.5;
	  }

	  /* Sum of the Fourier transformed hoppings */
	  cx_double ek1 = ts.ek1(kx, ky, kz);
	  
	  /* Fermi density */
	  double n_minus = 1.0;
	  double n_plus = 0.0;
	  if ( kB * T < 1e-15 ) {
	    n_minus = 1.0;
	    n_plus = 0.0;
	  } else {
	    double Em, Ep;
	    std::tie(Em, Ep) = calc_single_particle_energy2(ts, kx, ky, kz, delta);    
	    n_minus = fermi_density(Em, kB*T, mu);
	    n_plus = fermi_density(Ep, kB*T, mu);	    
	  }

	  sum += factor * (n_minus + n_plus);
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    // sum /= (double)( n_sites * delta );
    sum /= (double)( n_sites );    
  }

  return sum;
}

double calc_chemical_potential_bilayer3(path& base_dir, int L, hoppings_bilayer2 const& ts, double filling, double T, double delta, CubaParam const& cbp, bool continuous_k){  
  std::cout << "Finding the chemical potential for delta = " << delta << ", T = " << T << std::endl;
  double target = filling * 2;  // multiplied by 2
  double mu = 0.0;

  using std::placeholders::_1;
  auto eq = std::bind( elec_filling_eq_bilayer, L, std::ref(ts), _1, T, delta, std::ref(cbp), continuous_k );

  BinarySearch bs;
  // bool sol_found = bs.find_solution( mu, target, eq );

  // for check
  bool sol_found = bs.find_solution( mu, target, eq, false, 0.01, true );
  

  /* Output */
  std::cout << "mu = " << mu << std::endl;

  /* Precision */
  int prec = 12;
  
  ofstream out_mu;
  out_mu.open( base_dir / "chemical_potential.text");
  out_mu << "mu = " << std::setprecision(prec) << std::setw( prec ) << mu << std::endl;
  out_mu.close();

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

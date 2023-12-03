/*****************************************************************************
*
* Functions for finding the critical point.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "find_critical_point.h"
#include "BinarySearch.h"
#include "calc_chemical_potential.h"
#include "Hartree_Fock.h"

/* Instantiation */
FindQCPIntegrandBilayer fqcpib;

/* Integrand */
int find_critical_point_bilayer_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return fqcpib.calc(ndim, xx, ncomp, ff, userdata);
}

double self_consistent_eq_point_bilayer2(double tz, rpa::parameters& pr, CubaParam const& cbp){
  /* Setting the parameter */
  pr.t4 = tz;

  /* Extracting parameters */
  int L = pr.L;
  bool continuous_k = pr.continuous_k;

  /* Integral */
  double sum = 0;

  /* Hopping parameters */
  std::unique_ptr<hoppings_bilayer2> ts = hoppings_bilayer2::mk_bilayer3(pr);
  
  /* Changing the parameters */
  double mu = calc_chemical_potential_bilayer3(L, *ts, pr.filling, 0, 0, cbp, continuous_k, false);
  
  fqcpib.set_parameters(*ts, mu);
    
  if ( continuous_k ) {  
    /* For Cuba */
    double epsrel = 1e-10;
    double epsabs = 1e-10;
    int ncomp = 1;    
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

    /* Cuhre */
    Cuhre(cbp.NDIM, ncomp, find_critical_point_bilayer_integrand, cbp.userdata, cbp.nvec, epsrel, epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    // for(int comp = 0; comp < ncomp; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);      

    sum = integral[0] / 4.;
  } else { /* Integral for a finite-size system */
    double eps = 1e-12;
    double k1 = 2. * M_PI / (double)L;
  
    for(int z=0; z < 2; z++){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; x++){
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; y++){
	  double ky = k1 * y;      

	  /* Prefactor */
	  double factor = BZ_factor_square_half_filling(kx, ky);
	  if (factor == 0) continue;
	  
	  double qvec[3] = {kx, ky, kz};
	  sum += factor * fqcpib.integrand(qvec);
	} /* end for y */
      } /* end for x */
    } /* end for z */

    int n_sites = L * L * 2;  
    sum /= (double)(n_sites);    
  }

  return sum;
}

double find_critical_point_bilayer(rpa::parameters& pr){
  std::cout << "Finding a critical point for U = " << pr.U << std::endl;
  
  /* Getting parameters */
  double U = pr.U;
  assert(pr.filling == 0.5);  // Assume half filling  
  
  /* Parameters for Cuba */
  CubaParam cbp(pr);
  double target = 1. / U;

  /* Setting the parameter to an initial value. */
  double tz = 0;
  if ( pr.init_value == std::numeric_limits<double>::quiet_NaN() ) {
    tz = std::max(0.53*U, 1e-5);
  } else {
    tz = pr.init_value;
  }
  
  using std::placeholders::_1;
  auto scc = std::bind( self_consistent_eq_point_bilayer2, _1, std::ref(pr), std::ref(cbp) );

  BinarySearch bs(pr.continuous_k);
  bs.set_x_MIN(0);
  bool verbose = true;
  bool sol_found = bs.find_solution( tz, target, scc, false, 0.0, verbose );
  
  double QCP = 0;
  if ( sol_found ) {
    QCP = tz;
  } else {
    std::cerr << "Not found." << std::endl;
    QCP = 0;
  }
  
  return QCP;
}

void find_critical_point_bilayer_output(path& base_dir, rpa::parameters& pr){
  std::cout << "Finding the critical point..." << std::endl;
  double QCP = find_critical_point_bilayer(pr);
  
  /* Precision */
  int prec = 12;
  
  /* Output */
  std::cout << "QCP = " << QCP << std::endl;
  
  ofstream out_QCP;
  out_QCP.open( base_dir / "critical-point.text");
  out_QCP << "QCP = " << std::setprecision(prec) << std::setw( prec ) << QCP << std::endl;
  out_QCP.close();  
}

/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

// #include <functional>
#include <armadillo>
#include <cuba.h>
#include "calc_gap.h"
#include "rpa_util.h"
#include "calc_intensity.h"
#include "mat_elem.h"

/* Index: NCOMP = (ope * nelem + elem) | (real=0 or imag=1) */
std::tuple<int, int> comp_to_ope_elem(int comp, int nelem){
  comp >>= 1;
  int ope = comp / nelem;
  int elem = comp % nelem;
  return std::make_tuple(ope, elem);
}

/* Member functions of ResponseFuncIntegrand */
void ResponseFuncIntegrand::update_parameters(double _T, double _delta, double _mu, cx_double _omega, MatElemF const& _me_F){
  T_ = _T;
  delta_ = _delta;
  mu_ = _mu;
  omega_ = _omega;
  me_F_ = _me_F;  // without copying the tables
  me_F_.is_table_set_ = false;
};

/* Member functions of ResponseFuncIntegrandBilayer */
void ResponseFuncIntegrandBilayer::set_parameters(hoppings_bilayer2 const& h, double T, double delta, double mu, cx_double omega, MatElemF const& me_F){
  hb_ = h;
  ResponseFuncIntegrand::update_parameters(T, delta, mu, omega, me_F);
}

int ResponseFuncIntegrandBilayer::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) const {
   /* Reset */
   for(int comp=0; comp<*ncomp; comp++){ ff[comp] = 0; }
   				    
   /* Wavenumbers */
   double k1 = xx[0] * 2 * M_PI;
   double k2 = xx[1] * 2 * M_PI;
   
   double kx = 0.5 * (k2 + k1);
   double ky = 0.5 * (k2 - k1);

   /* Polarizaiton */
   cx_double F00k[NSUBL*NSUBL];   
   cx_double Fpmk[NSUBL*NSUBL];
   cx_double Fzzk[NSUBL*NSUBL];
       
   /* Sum over kz */
   for(int z=0; z < 2; z++){       
     double kz = M_PI * z;	  
     
     /* Sum over bands (signs) */
     for(int sg1=-1; sg1<=1; sg1+=2){
       for(int sg2=-1; sg2<=1; sg2+=2){
	 // int sg2 = - sg1; /* Opposite sign */
	 cx_double prefactor = calc_prefactor_bare_res_func_bilayer(sg1, sg2, *ts(), T(), kx, ky, kz, me_F()->qx(), me_F()->qy(), me_F()->qz(), omega(), delta(), mu());
    
	 /* Getting the polarization */
	 me_F()->calc_mat_elems(hb_, delta(), kx, ky, kz, sg1, sg2, F00k, Fpmk, Fzzk);

	 /* Integrand */
	 for(int comp=0; comp<*ncomp; comp+=2){
	   int ope, elem;
	   std::tie(ope, elem) = comp_to_ope_elem(comp, NSUBL*NSUBL);
	   if ( ope == 0 ) {  /* Fpmk */
	     ff[comp] += std::real(prefactor * Fpmk[elem]);
	     ff[comp+1] += std::imag(prefactor * Fpmk[elem]);
	   } else /* ope == 1 */ { /* Fzzk */	      
	     ff[comp] += std::real(prefactor * Fzzk[elem]);
	     ff[comp+1] += std::imag(prefactor * Fzzk[elem]);
	   }
	 }
       }
     }
   }
   
   return 0;   
}


cx_double calc_intensity_bilayer(int L, hoppings const& ts, double mu, double U, double delta, double qx, double qy, double qz, cx_double omega, bool zz){
  double k1 = 2. * M_PI / (double)L;
  
  cx_double A = 0, B = 0, C = 0, D = 0;
  for(int z=0; z < 2; z++){    
    double kz = M_PI * z;
    for(int x=-L/2; x < L/2; x++){    
      double kx = k1 * x;
      for(int y=-L/2; y < L/2; y++){
	double ky = k1 * y;
	add_to_sus_mat2( ts, mu, A, B, C, D, qx, qy, qz, kx, ky, kz, delta, omega, zz );	
      }
    }
  }

  int n_sites = L * L * 2;
  double norm = 2. / (double)n_sites;
  A *= norm;
  B *= norm;
  C *= norm;
  D *= norm;
  
  /* RPA */
  arma::cx_mat chi0_mat(2,2);
  chi0_mat(0,0) = A;   // (A, A) correlation
  chi0_mat(0,1) = B;   // (A, B)
  chi0_mat(1,0) = C;   // (B, A)
  chi0_mat(1,1) = D;   // (B, B)

  /* Transverse = < \sigma^- \sigma^+ >; Longitudinal (zz) = < \sigma^z \sigma^z > */
  /* Note that 2 < \sigma^- \sigma^+ > = < \sigma^z \sigma^z > (U -> 0 for the SU(2) case) */  
  double factor_channel = 1.0;
  if ( zz ) { factor_channel = 0.5; }
  
  arma::cx_mat denom = arma::eye<arma::cx_mat>(NSUBL,NSUBL) - factor_channel * U * chi0_mat;
  arma::cx_mat chi_mat = chi0_mat * arma::inv(denom);
  
  // sigma-to-spin factor
  double factor_operator = 0.5;  
  cx_double chi = factor_operator * factor_operator * ( chi_mat(0,0) - chi_mat(1,0) - chi_mat(0,1) + chi_mat(1,1) );  

  // // for check
  // std::cout << chi0_mat << std::endl;
  // std::cout << chi_mat << std::endl;
  
  return chi;  
}

/* Instantiation */
ResponseFuncIntegrandBilayer rfib;

int integrand_wrapper(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return rfib.calc(ndim, xx, ncomp, ff, userdata);
}

std::tuple<arma::cx_mat, arma::cx_mat> calc_bare_response_bilayer(int L, hoppings_bilayer2 const& ts, double mu, double U, double T, double delta, CubaParam const& cbp, MatElemF const& me_F, cx_double omega, bool continuous_k){
  arma::cx_mat chi0_pm(NSUBL,NSUBL,arma::fill::zeros);
  arma::cx_mat chi0_zz_u(NSUBL,NSUBL,arma::fill::zeros);
  
  if ( continuous_k ) {
    /* Changing the parameters */
    rfib.set_parameters(ts, T, delta, mu, omega, me_F);

    /* For Cuba */
    int nregions, neval, fail;
    cubareal integral[cbp.NCOMP], error[cbp.NCOMP], prob[cbp.NCOMP];
    
    /* Cuhre */
    Cuhre(cbp.NDIM, cbp.NCOMP, integrand_wrapper, cbp.userdata, cbp.nvec, cbp.epsrel, cbp.epsabs, cbp.flags, cbp.mineval, cbp.maxeval, cbp.key, cbp.statefile, cbp.spin, &nregions, &neval, &fail, integral, error, prob);
    
    // for check
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    // for(int comp = 0; comp < cbp.NCOMP; comp++ )
    //   printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    // 	     (double)integral[comp], (double)error[comp], (double)prob[comp]);      
      
    /* Extracting the integral results */
    for(int comp=0; comp<cbp.NCOMP; comp+=2){
      int ope, elem;
      std::tie(ope, elem) = comp_to_ope_elem(comp, NSUBL*NSUBL);
      int g1 = elem / NSUBL;
      int g2 = elem % NSUBL;
      cx_double int_ope_elem(integral[comp], integral[comp+1]);
      if ( ope == 0 ) {  /* Ppmk */
	chi0_pm(g1,g2) = int_ope_elem;
      } else /* ope == 1 */ { /* me_Fzk */	      
	chi0_zz_u(g1,g2) = int_ope_elem;
      }
    }

    chi0_pm /= 2; /* A factor (2*M_PI) cancels because of the scale change of the integration variables. */
    chi0_zz_u /= 2;
  } else { /* Integral for a finite-size system */
    double k1 = 2. * M_PI / L;
    for(int z=0; z < 2; z++){    
      double kz = M_PI * z;
      for(int x=-L/2; x < L/2; x++){    
	double kx = k1 * x;
	for(int y=-L/2; y < L/2; y++){
	  double ky = k1 * y;
	  add_to_sus_mat4( ts, T, mu, chi0_pm, chi0_zz_u, kx, ky, kz, me_F, delta, omega );	  
	}
      }
    }
    
    int n_units = L * L;  // Number of unit cells
    chi0_pm /= (double)(n_units);
    chi0_zz_u /= (double)(n_units);
  }

  /* Adding the contribution from the down spin */
  arma::cx_mat chi0_zz_d(chi0_zz_u);
  /* Assume NSUBL == 2. */  
  if ( NSUBL == 2 ) {
    chi0_zz_d.swap_rows( 0, 1 );
    chi0_zz_d.swap_cols( 0, 1 );
  } else {
    std::cerr << "NSUBL is not 2.\n";
    std::exit(EXIT_FAILURE);
  }
  arma::cx_mat chi0_zz = chi0_zz_u + chi0_zz_d;
  return std::make_tuple(chi0_pm, chi0_zz);
}

std::tuple<cx_double, cx_double> calc_intensity_bilayer2(int L, hoppings_bilayer2& ts, double mu, double U, double T, double delta, CubaParam const& cbp, MatElemF const& me_F, cx_double omega, bool continuous_k){
  arma::cx_mat chi0_pm(NSUBL,NSUBL,arma::fill::zeros);
  arma::cx_mat chi0_zz(NSUBL,NSUBL,arma::fill::zeros);

  /* Calculating the bare response functions */
  std::tie(chi0_pm, chi0_zz) = calc_bare_response_bilayer(L, ts, mu, U, T, delta, cbp, me_F, omega, continuous_k);
  
  /* RPA */
  /* Transverse = < \sigma^- \sigma^+ >; Longitudinal (zz) = < \sigma^z \sigma^z > */
  /* Note that 2 < \sigma^- \sigma^+ > = < \sigma^z \sigma^z > (U -> 0 for the SU(2) case) */  
  arma::cx_mat denom_pm = arma::eye<arma::cx_mat>(NSUBL,NSUBL) - U * chi0_pm;
  arma::cx_mat denom_zz = arma::eye<arma::cx_mat>(NSUBL,NSUBL) - 0.5 * U * chi0_zz;  
  arma::cx_mat chi_pm = arma::inv(denom_pm) * chi0_pm;
  arma::cx_mat chi_zz = arma::inv(denom_zz) * chi0_zz;
  
  // sigma-to-spin factor
  cx_double chi_xy = 0.5 * arma::accu(chi_pm);
  cx_double chi_z = 0.25 * arma::accu(chi_zz);  

  // // for check
  // std::cout << chi0_pm << std::endl;
  // std::cout << chi_pm << std::endl;  
  // std::cout << chi0_zz << std::endl;
  // std::cout << chi_zz << std::endl;
  
  return std::make_tuple(chi_xy, chi_z);
}

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

/* Index: NCOMP = (ope * nelem + elem) | (real=0 or imag=1) */
std::tuple<int, int> comp_to_ope_elem(int comp, int nelem){
  comp >>= 1;
  int ope = comp / nelem;
  int elem = comp % nelem;
  return std::make_tuple(ope, elem);
}

class Integrand_k {
public:
  Integrand_k():ts_(),Pz_(){}
  void set_parameters(hoppings_bilayer2 const& ts, double delta, cx_double omega, Polarization const& Pz_pt);
  int calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);
  double delta() const { return delta_; }
  cx_double omega() const { return omega_; }
  
private:  
  hoppings_bilayer2 ts_;
  double delta_;
  cx_double omega_;
  Polarization Pz_;
};

void Integrand_k::set_parameters(hoppings_bilayer2 const& ts, double delta, cx_double omega, Polarization const& Pz){
  ts_ = ts;
  delta_ = delta;
  omega_ = omega;
  Pz_ = Pz;
}

int Integrand_k::calc(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
   /* Reset */
   for(int comp=0; comp<*ncomp; comp++){ ff[comp] = 0; }
   				    
   /* Wavenumbers */
   double k1 = xx[0] * 2 * M_PI;
   double k2 = xx[1] * 2 * M_PI;
   
   double kx = 0.5 * (k2 + k1);
   double ky = 0.5 * (k2 - k1);
   double kx2 = kx + Pz_.qx();
   double ky2 = ky + Pz_.qy();

   /* Polarizaiton */
   cx_double Ppmk[NSUBL*NSUBL];
   cx_double Pzzk[NSUBL*NSUBL];
       
   /* Sum over kz */
   for(int z=0; z < 2; z++){       
     double kz = M_PI * z;	  
     double kz2 = kz + Pz_.qz();
   
     /* Energies */
     cx_double ek1 = ts_.ek1(kx, ky, kz);
     cx_double ek23 = ts_.ek23(kx, ky, kz);
     cx_double ekz = ts_.ekz(kx, ky, kz);
     cx_double ek_q1 = ts_.ek1( kx2, ky2, kz2 );
     cx_double ek_q23 = ts_.ek23( kx2, ky2, kz2 );
     cx_double ek_qz = ts_.ekz( kx2, ky2, kz2 );  

     /* Sum over bands (signs) */
     for(int sg1=-1; sg1<=1; sg1+=2){
       int sg2 = - sg1; /* Opposite sign */
       int sg1i = (sg1+1) >> 1;    
       int sg2i = (sg2+1) >> 1;
   
       /* Electron density at zero temperature */
       int n1 = ( ( sg1 + 1 ) >> 1 ) ^ 1;  /* -1 -> 1, 1 -> 0 */
       int n2 = n1 ^ 1;
       double n_diff = (double)(n2 - n1);

       /* Prefactor */
       double Ek = eigenenergy_HF(sg1, ek1, ek23, ekz, ts_.tz, kz, delta());
       double Ek_q = eigenenergy_HF(sg2, ek_q1, ek_q23, ek_qz, ts_.tz, kz2, delta());
       cx_double diff_E = omega() - (Ek_q - Ek);    
       cx_double prefactor = n_diff / diff_E;
       
       /* Getting the polarization */
       Pz_.calc_polarization(ts_, delta(), kx, ky, kz, sg1, sg2, Ppmk, Pzzk);

       /* Integrand */
       for(int comp=0; comp<*ncomp; comp+=2){
	 int ope, elem;
	 std::tie(ope, elem) = comp_to_ope_elem(comp, NSUBL*NSUBL);
	 if ( ope == 0 ) {  /* Ppmk */
	   ff[comp] += std::real(prefactor * Ppmk[elem]);
	   ff[comp+1] += std::imag(prefactor * Ppmk[elem]);
	 } else /* ope == 1 */ { /* Pzzk */	      
	   ff[comp] += std::real(prefactor * Pzzk[elem]);
	   ff[comp+1] += std::imag(prefactor * Pzzk[elem]);
	 }
       }
     }
   }
   
   return 0;   
}

/* Creating an instance */
Integrand_k integrand_k;

int integrand_wrapper(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata){
  return integrand_k.calc(ndim, xx, ncomp, ff, userdata);
}

std::tuple<cx_double, cx_double> calc_intensity_bilayer2(int L, hoppings_bilayer2& ts, double mu, double U, double delta, CubaParam const& cbp, double qx, double qy, double qz, Polarization& Pz, cx_double omega, bool continuous_k){
  arma::cx_mat chi0_pm(NSUBL,NSUBL,arma::fill::zeros);
  arma::cx_mat chi0_zz_u(NSUBL,NSUBL,arma::fill::zeros);      
  
  if ( continuous_k ) {
    /* Changing the parameters */
    integrand_k.set_parameters(ts, delta, omega, Pz);

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
      } else /* ope == 1 */ { /* Pzzk */	      
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
	  add_to_sus_mat4( ts, mu, chi0_pm, chi0_zz_u, qx, qy, qz, kx, ky, kz, Pz, delta, omega );
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
  
  /* RPA */
  /* Transverse = < \sigma^- \sigma^+ >; Longitudinal (zz) = < \sigma^z \sigma^z > */
  /* Note that 2 < \sigma^- \sigma^+ > = < \sigma^z \sigma^z > (U -> 0 for the SU(2) case) */  
  arma::cx_mat denom_pm = arma::eye<arma::cx_mat>(NSUBL,NSUBL) - U * chi0_pm;
  arma::cx_mat denom_zz = arma::eye<arma::cx_mat>(NSUBL,NSUBL) - 0.5 * U * chi0_zz;  
  arma::cx_mat chi_pm = chi0_pm * arma::inv(denom_pm);
  arma::cx_mat chi_zz = chi0_zz * arma::inv(denom_zz);  
  
  // sigma-to-spin factor
  cx_double chi_xy = 0.25 * arma::accu(chi_pm);
  cx_double chi_z = 0.25 * arma::accu(chi_zz);  

  // for check
  std::cout << chi0_pm << std::endl;
  std::cout << chi_pm << std::endl;  
  std::cout << chi0_zz << std::endl;
  std::cout << chi_zz << std::endl;
  
  return std::make_tuple(chi_xy, chi_z);
}

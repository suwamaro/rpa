/*****************************************************************************
*
* Functions for calculating the Raman spectrum
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "rpa_util.h"
#include "Hartree_Fock.h"
#include "mat_elem.h"

/* Member functions of MatElem */
MatElem::MatElem(int Lx, int Ly, int Lz, int nsub){
  is_table_set_ = false;
  Lx_ = Lx;
  Ly_ = Ly;
  Lz_ = Lz;
  nsub_ = nsub;
}

void MatElem::set_q(double qx, double qy, double qz){
  qx_ = qx;
  qy_ = qy;
  qz_ = qz;
}

/* Assume that k = 2pi / L * m, where m is an integer. */
int MatElem::pullback(double k, int L) const {
  int l = rint(k * L / (2.*M_PI));
  while ( l < 0 ) l += L;
  while ( l >= L ) l -= L;
  /* Return 0 <= l < L */
  return l;
}
  
int MatElem::xyz_to_index(int x, int y, int z) const {
  return z * Lx() * Ly() + y * Lx() + x;
}
  
int MatElem::k_to_index(double kx, double ky, double kz) const {    
  int lx = pullback(kx, Lx());
  int ly = pullback(ky, Ly());
  int lz = pullback(kz, Lz());    
  return xyz_to_index(lx, ly, lz);
}


/* Member functions of MatElemF */
MatElemF::MatElemF(int Lx, int Ly, int Lz, int nsub):MatElem(Lx, Ly, Lz, nsub){}

void MatElemF::set_table(hoppings2 const& ts, double delta){
  F00_ = new cx_double[table_size()];  
  Fpm_ = new cx_double[table_size()];
  Fzz_ = new cx_double[table_size()];  
  build_table(ts, delta);
  is_table_set_ = true;
}

std::size_t MatElemF::table_size() const { return Lx()*Ly()*Lz()*nbands()*nbands()*nsub()*nsub(); }
  
void MatElemF::calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* F00, cx_double* Fpm, cx_double* Fzz) const {
  /* Adding the results */
  cx_double ek1 = ts.ek1(kx, ky, kz);
  double kx2 = kx + qx();
  double ky2 = ky + qy();
  double kz2 = kz + qz();
  cx_double ek_q1 = ts.ek1( kx2, ky2, kz2 );
  cx_double tz = ts.tz;

  cx_vec X1up = gs_HF1(up_spin, sg1, ek1, tz, kz, delta);
  cx_vec X1down = gs_HF1(down_spin, sg1, ek1, tz, kz, delta);
  cx_vec X2up = gs_HF1(up_spin, sg2, ek_q1, tz, kz2, delta);
  cx_vec X2down = gs_HF1(down_spin, sg2, ek_q1, tz, kz2, delta);

  cx_mat P1up = arma::kron(arma::trans(arma::cx_mat(X1up)), X1up);
  // cx_mat P1down = arma::kron(arma::trans(arma::cx_mat(X1down)), X1down);
  cx_mat P2up = arma::kron(arma::trans(arma::cx_mat(X2up)), X2up);
  cx_mat P2down = arma::kron(arma::trans(arma::cx_mat(X2down)), X2down);    
  
  mat tau_A = mat({1., 0., 0., 0.});
  mat tau_B = mat({0., 0., 0., 1.});
  tau_A.set_size(2,2);  
  tau_B.set_size(2,2);
  
  cx_mat Pauli_0 = cx_mat({1., 0., 0., 1.});
  cx_mat Pauli_p = cx_mat({0., 0., 1., 0.});
  cx_mat Pauli_m = cx_mat({0., 1., 0., 0.});
  cx_mat Pauli_z = cx_mat({1., 0., 0., - 1.});
  Pauli_0.set_size(2,2);
  Pauli_p.set_size(2,2);
  Pauli_m.set_size(2,2);  
  Pauli_z.set_size(2,2);
  
  cx_mat Sigma_A0 = arma::kron(tau_A, Pauli_0);
  cx_mat Sigma_Ap = arma::kron(tau_A, Pauli_p);
  cx_mat Sigma_Am = arma::kron(tau_A, Pauli_m);
  cx_mat Sigma_Az = arma::kron(tau_A, Pauli_z);
  cx_mat Sigma_B0 = arma::kron(tau_B, Pauli_0);
  cx_mat Sigma_Bp = arma::kron(tau_B, Pauli_p);
  cx_mat Sigma_Bm = arma::kron(tau_B, Pauli_m);
  cx_mat Sigma_Bz = arma::kron(tau_B, Pauli_z);      

  cx_double F_A0_A0 = arma::trace(P1up * Sigma_A0 * P2up * Sigma_A0);
  cx_double F_A0_B0 = arma::trace(P1up * Sigma_A0 * P2up * Sigma_B0);
  cx_double F_B0_A0 = arma::trace(P1up * Sigma_B0 * P2up * Sigma_A0);
  cx_double F_B0_B0 = arma::trace(P1up * Sigma_B0 * P2up * Sigma_B0);    
  F00[0] = F_A0_A0;
  F00[1] = F_A0_B0;
  F00[2] = F_B0_A0;
  F00[3] = F_B0_B0;
  
  cx_double F_Az_Az = arma::trace(P1up * Sigma_Az * P2up * Sigma_Az);
  cx_double F_Az_Bz = arma::trace(P1up * Sigma_Az * P2up * Sigma_Bz);
  cx_double F_Bz_Az = arma::trace(P1up * Sigma_Bz * P2up * Sigma_Az);
  cx_double F_Bz_Bz = arma::trace(P1up * Sigma_Bz * P2up * Sigma_Bz);    
  Fzz[0] = F_Az_Az;
  Fzz[1] = F_Az_Bz;
  Fzz[2] = F_Bz_Az;
  Fzz[3] = F_Bz_Bz;
  
  cx_double F_Ap_Am = arma::trace(P1up * Sigma_Ap * P2down * Sigma_Am);
  cx_double F_Ap_Bm = arma::trace(P1up * Sigma_Ap * P2down * Sigma_Bm);
  cx_double F_Bp_Am = arma::trace(P1up * Sigma_Bp * P2down * Sigma_Am);
  cx_double F_Bp_Bm = arma::trace(P1up * Sigma_Bp * P2down * Sigma_Bm);  
  Fpm[0] = F_Ap_Am;
  Fpm[1] = F_Ap_Bm;
  Fpm[2] = F_Bp_Am;
  Fpm[3] = F_Bp_Bm;
  
  // cx_vec X1up = gs_HF1(up_spin, sg1, ek1, tz, kz, delta);
  // // cx_vec X1down = gs_HF1(down_spin, sg1, ek1, tz, kz, delta);
  // cx_vec X2up = gs_HF1(up_spin, sg2, ek_q1, tz, kz2, delta);  
  // cx_vec X2down = gs_HF1(down_spin, sg2, ek_q1, tz, kz2, delta);
  
  // cx_double F_zu_g1 = arma::cdot(X1up, opek.Gamma_1z * X2up);
  // cx_double F_zu_gm1 = arma::cdot(X1up, opek.Gamma_2z * X2up);
  	    
  // /* Assume the number of sublattices is 2. */
  // pzz[ 0 ] = std::norm(F_zu_g1);      
  // pzz[ 1 ] = F_zu_g1 * std::conj(F_zu_gm1);
  // pzz[ 2 ] = F_zu_gm1 * std::conj(F_zu_g1);
  // pzz[ 3 ] = std::norm(F_zu_gm1);      

  /* <sigma_minus sigma_plus> */  
  // cx_double F_m_g1 = arma::cdot(X1down, opek.Gamma_1minus * X2up);  
  // cx_double F_m_gm1 = arma::cdot(X1down, opek.Gamma_2minus * X2up);
  // ppm[ 0 ] = std::norm(F_m_g1);  
  // ppm[ 1 ] = F_m_g1 * std::conj(F_m_gm1);
  // ppm[ 2 ] = std::conj(ppm[1]);
  // ppm[ 3 ] = std::norm(F_m_gm1);

  // cx_double F_p_g1 = arma::cdot(X1up, opek.Gamma_1plus * X2down);
  // cx_double F_p_gm1 = arma::cdot(X1up, opek.Gamma_2plus * X2down);    
  // ppm[ 0 ] = std::norm(F_p_g1);      
  // ppm[ 1 ] = F_p_g1 * std::conj(F_p_gm1);
  // ppm[ 2 ] = F_p_gm1 * std::conj(F_p_g1);  
  // ppm[ 3 ] = std::norm(F_p_gm1);
}

void MatElemF::build_table(hoppings2 const& ts, double delta){
  cx_double F00[nsub()*nsub()];  
  cx_double Fpm[nsub()*nsub()];
  cx_double Fzz[nsub()*nsub()];
  for(int x=0; x < Lx(); x++){    
    double kx = 2. * M_PI / Lx() * x;
    for(int y=0; y < Ly(); y++){    
      double ky = 2. * M_PI / Ly() * y;
      for(int z=0; z < Lz(); z++){    
	double kz = 2. * M_PI / Lz() * z;
	std::size_t xyz_idx = xyz_to_index(x,y,z) * nbands() * nbands() * nsub() * nsub();
	
	for(int sg1=-1; sg1<=1; sg1+=2){
	  for(int sg2=-1; sg2<=1; sg2+=2){
	    // int sg2 = - sg1; /* Opposite sign */
	    calc_mat_elems(ts, delta, kx, ky, kz, sg1, sg2, F00, Fpm, Fzz);
	  
	    int sg1i = (sg1+1) >> 1;
	    int sg2i = (sg2+1) >> 1;
	    std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
	    std::size_t xyz_bands_idx = xyz_idx + bands_idx;
	    memcpy(F00_ + xyz_bands_idx, F00, sizeof(cx_double)*nsub()*nsub() );	    
	    memcpy(Fpm_ + xyz_bands_idx, Fpm, sizeof(cx_double)*nsub()*nsub() );
	    memcpy(Fzz_ + xyz_bands_idx, Fzz, sizeof(cx_double)*nsub()*nsub() );
	  }
	}
      }
    }
  }
}

void MatElemF::get_00(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const {
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( F00k, F00_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}

void MatElemF::get_pm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fpmk) const {
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( Fpmk, Fpm_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}

void MatElemF::get_zz(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const {    
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( Fzzk, Fzz_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}
  
MatElemF::~MatElemF(){
  if ( F00_ != nullptr ) {
    delete[] F00_;
  }
  
  if ( Fpm_ != nullptr ) {
    delete[] Fpm_;
  }

  if ( Fzz_ != nullptr ) {
    delete[] Fzz_;
  }
}

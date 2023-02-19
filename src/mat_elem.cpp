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
#include "calc_Raman.h"

/* Member functions of MatElem */
MatElem::MatElem(int Lx, int Ly, int Lz){
  is_table_set_ = false;
  Lx_ = Lx;
  Ly_ = Ly;
  Lz_ = Lz;
}

MatElem::MatElem(int Lx, int Ly, int Lz, int n_coefs){
  is_table_set_ = false;
  Lx_ = Lx;
  Ly_ = Ly;
  Lz_ = Lz;
  n_coefs_ = n_coefs;
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


/* Member functions of BasicMatrix */
BasicMatrix::BasicMatrix(){
  tau_A = mat({1., 0., 0., 0.});
  tau_B = mat({0., 0., 0., 1.});
  tau_A.set_size(2,2);  
  tau_B.set_size(2,2);
  
  Pauli_0 = cx_mat({1., 0., 0., 1.});
  Pauli_p = cx_mat({0., 0., 1., 0.});
  Pauli_m = cx_mat({0., 1., 0., 0.});
  Pauli_z = cx_mat({1., 0., 0., - 1.});
  Pauli_0.set_size(2,2);
  Pauli_p.set_size(2,2);
  Pauli_m.set_size(2,2);  
  Pauli_z.set_size(2,2);
  
  Sigma_A0 = arma::kron(tau_A, Pauli_0);
  Sigma_Ap = arma::kron(tau_A, Pauli_p);
  Sigma_Am = arma::kron(tau_A, Pauli_m);
  Sigma_Az = arma::kron(tau_A, Pauli_z);
  Sigma_B0 = arma::kron(tau_B, Pauli_0);
  Sigma_Bp = arma::kron(tau_B, Pauli_p);
  Sigma_Bm = arma::kron(tau_B, Pauli_m);
  Sigma_Bz = arma::kron(tau_B, Pauli_z);
}


/* Member functions of MatElemF */
MatElemF::MatElemF(int Lx, int Ly, int Lz, int nsub):MatElem(Lx, Ly, Lz),BasicMatrix(),nsub_(nsub){}

void MatElemF::set_table(hoppings2 const& ts, double delta){
  if ( F00_up_ != nullptr ) {
    delete[] F00_up_;
  }
  F00_up_ = new cx_double[table_size()];

  if ( F00_down_ != nullptr ) {
    delete[] F00_down_;
  }
  F00_down_ = new cx_double[table_size()];  

  if ( Fpm_ != nullptr ) {
    delete[] Fpm_;
  }  
  Fpm_ = new cx_double[table_size()];

  if ( Fzz_up_ != nullptr ) {
    delete[] Fzz_up_;
  }  
  Fzz_up_ = new cx_double[table_size()];

  if ( Fzz_down_ != nullptr ) {
    delete[] Fzz_down_;
  }  
  Fzz_down_ = new cx_double[table_size()];  
  
  build_table(ts, delta);
  assign_is_table_set(true);
}

std::size_t MatElemF::table_size() const { return Lx()*Ly()*Lz()*nbands()*nbands()*nsub()*nsub(); }
  
void MatElemF::calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* F00_up, cx_double* F00_down, cx_double* Fpm, cx_double* Fzz_up, cx_double* Fzz_down) const {
  /* Adding the results */
  cx_double ek1 = ts.ek1(kx, ky, kz);
  double kx2 = kx + qx();
  double ky2 = ky + qy();
  double kz2 = kz + qz();
  cx_double ek_q1 = ts.ek1(kx2, ky2, kz2);
  cx_double tz = ts.tz;

  cx_vec X1up = gs_HF1(up_spin, sg1, ek1, tz, kz, delta);
  cx_vec X1down = gs_HF1(down_spin, sg1, ek1, tz, kz, delta);
  cx_vec X2up = gs_HF1(up_spin, sg2, ek_q1, tz, kz2, delta);
  cx_vec X2down = gs_HF1(down_spin, sg2, ek_q1, tz, kz2, delta);

  cx_double A0_12up = arma::cdot(X1up, Sigma_A0 * X2up);  
  cx_double A0_21up = arma::cdot(X2up, Sigma_A0 * X1up);
  cx_double B0_12up = arma::cdot(X1up, Sigma_B0 * X2up);
  cx_double B0_21up = arma::cdot(X2up, Sigma_B0 * X1up);  
  F00_up[0] = A0_21up * A0_12up;
  F00_up[1] = A0_21up * B0_12up;
  F00_up[2] = B0_21up * A0_12up;
  F00_up[3] = B0_21up * B0_12up;

  cx_double A0_12down = arma::cdot(X1down, Sigma_A0 * X2down);  
  cx_double A0_21down = arma::cdot(X2down, Sigma_A0 * X1down);
  cx_double B0_12down = arma::cdot(X1down, Sigma_B0 * X2down);
  cx_double B0_21down = arma::cdot(X2down, Sigma_B0 * X1down);  
  F00_down[0] = A0_21down * A0_12down;
  F00_down[1] = A0_21down * B0_12down;
  F00_down[2] = B0_21down * A0_12down;
  F00_down[3] = B0_21down * B0_12down;

  cx_double Az_12up = arma::cdot(X1up, Sigma_Az * X2up);    
  cx_double Az_21up = arma::cdot(X2up, Sigma_Az * X1up);
  cx_double Bz_12up = arma::cdot(X1up, Sigma_Bz * X2up);  
  cx_double Bz_21up = arma::cdot(X2up, Sigma_Bz * X1up);
  Fzz_up[0] = Az_21up * Az_12up;
  Fzz_up[1] = Az_21up * Bz_12up;
  Fzz_up[2] = Bz_21up * Az_12up;
  Fzz_up[3] = Bz_21up * Bz_12up;

  cx_double Az_12down = arma::cdot(X1down, Sigma_Az * X2down);    
  cx_double Az_21down = arma::cdot(X2down, Sigma_Az * X1down);
  cx_double Bz_12down = arma::cdot(X1down, Sigma_Bz * X2down);  
  cx_double Bz_21down = arma::cdot(X2down, Sigma_Bz * X1down);
  Fzz_down[0] = Az_21down * Az_12down;
  Fzz_down[1] = Az_21down * Bz_12down;
  Fzz_down[2] = Bz_21down * Az_12down;
  Fzz_down[3] = Bz_21down * Bz_12down;

  cx_double Ap_21 = arma::cdot(X2up, Sigma_Ap * X1down);
  cx_double Bp_21 = arma::cdot(X2up, Sigma_Bp * X1down);
  Fpm[0] = std::norm(Ap_21);
  Fpm[1] = Ap_21 * std::conj(Bp_21);
  Fpm[2] = Bp_21 * std::conj(Ap_21);
  Fpm[3] = std::norm(Bp_21);
}

void MatElemF::build_table(hoppings2 const& ts, double delta){
  cx_double F00_up[nsub()*nsub()];
  cx_double F00_down[nsub()*nsub()];    
  cx_double Fpm[nsub()*nsub()];
  cx_double Fzz_up[nsub()*nsub()];
  cx_double Fzz_down[nsub()*nsub()];  
  for(int z=0; z < Lz(); ++z){    
    double kz = 2. * M_PI / Lz() * z;
    for(int y=0; y < Ly(); ++y){    
      double ky = 2. * M_PI / Ly() * y;	
      for(int x=0; x < Lx(); ++x){    
	double kx = 2. * M_PI / Lx() * x;
	std::size_t xyz_idx = xyz_to_index(x,y,z) * nbands() * nbands() * nsub() * nsub();
	
	for(int sg1=-1; sg1<=1; sg1+=2){
	  for(int sg2=-1; sg2<=1; sg2+=2){
	    calc_mat_elems(ts, delta, kx, ky, kz, sg1, sg2, F00_up, F00_down, Fpm, Fzz_up, Fzz_down);
	  
	    int sg1i = (sg1+1) >> 1;
	    int sg2i = (sg2+1) >> 1;
	    std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
	    std::size_t xyz_bands_idx = xyz_idx + bands_idx;
	    memcpy(F00_up_ + xyz_bands_idx, F00_up, sizeof(cx_double)*nsub()*nsub() );
	    memcpy(F00_down_ + xyz_bands_idx, F00_down, sizeof(cx_double)*nsub()*nsub() );	    	    
	    memcpy(Fpm_ + xyz_bands_idx, Fpm, sizeof(cx_double)*nsub()*nsub() );
	    memcpy(Fzz_up_ + xyz_bands_idx, Fzz_up, sizeof(cx_double)*nsub()*nsub() );
	    memcpy(Fzz_down_ + xyz_bands_idx, Fzz_down, sizeof(cx_double)*nsub()*nsub() );	    
	  }
	}
      }
    }
  }
}

void MatElemF::get_00_up(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const {
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( F00k, F00_up_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}

void MatElemF::get_00_down(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* F00k) const {
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( F00k, F00_down_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}

void MatElemF::get_pm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fpmk) const {
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( Fpmk, Fpm_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}

void MatElemF::get_zz_up(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const {    
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( Fzzk, Fzz_up_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}

void MatElemF::get_zz_down(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Fzzk) const {    
  std::size_t xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  std::size_t bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  std::size_t xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( Fzzk, Fzz_down_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}
  
MatElemF::~MatElemF(){
  if ( F00_up_ != nullptr ) {
    delete[] F00_up_;
  }

  if ( F00_down_ != nullptr ) {
    delete[] F00_down_;
  }  
  
  if ( Fpm_ != nullptr ) {
    delete[] Fpm_;
  }

  if ( Fzz_up_ != nullptr ) {
    delete[] Fzz_up_;
  }

  if ( Fzz_down_ != nullptr ) {
    delete[] Fzz_down_;
  }  
}


/* Member functions of MatElemK */
MatElemK::MatElemK(int Lx, int Ly, int Lz, int mu, int nu, std::vector<BondDelta> const& bonds):MatElem(Lx, Ly, Lz, 2),mu_(mu),nu_(nu),bonds_(bonds){}  // A factor of 2 comes from the two cases of the perturbation.

void MatElemK::set_table(hoppings2 const& ts, double delta){
  if ( R_up_ != nullptr ) {
    delete[] R_up_;
  }
  R_up_ = new cx_double[table_size()];
  
  build_table(ts, delta);
  assign_is_table_set(true);
}

std::size_t MatElemK::table_size() const { return Lx()*Ly()*Lz()*n_coefs(); }

void MatElemK::calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *R) const {
  /* The photon momenta are approximately set to zero. */
  vec3 ki {0, 0, 0}, kf {0, 0, 0};
  
  /* Formulation using up spin */
  int spin = up_spin;
  
  /* Two cases of the perturbation */
  cx_double R1 = calc_coef_eff_Raman_resonant(Lx(), ts, delta, kx, ky, kz, spin, bonds_, mu_, nu_, ki, kf);
  cx_double R2 = calc_coef_eff_Raman_resonant(Lx(), ts, delta, kx, ky, kz, spin, bonds_, nu_, mu_, - kf, - ki);
  
  R[0] = R1;
  R[1] = R2;

  // // for check
  // if ( mu_.z == 1 && nu_.z == 1 ) {
  //   cx_double R1_up = calc_coef_eff_Raman_resonant(Lx(), ts, delta, kx, ky, kz, up_spin, bonds_, mu_, nu_, ki, kf);
  //   cx_double R2_up = calc_coef_eff_Raman_resonant(Lx(), ts, delta, kx, ky, kz, up_spin, bonds_, nu_, mu_, - kf, - ki);
  //   cx_double R1_down = calc_coef_eff_Raman_resonant(Lx(), ts, delta, kx, ky, kz, down_spin, bonds_, mu_, nu_, ki, kf);
  //   cx_double R2_down = calc_coef_eff_Raman_resonant(Lx(), ts, delta, kx, ky, kz, down_spin, bonds_, nu_, mu_, - kf, - ki);

  //   std::cout << kx << "  " << ky << "  " << kz << "     " << R1_up << "   " << R1_down << "     " << R2_up << "   " << R2_down << std::endl;
  // }
}

void MatElemK::build_table(hoppings2 const& ts, double delta){
  for(int z=0; z < Lz(); ++z){    
    double kz = 2. * M_PI / Lz() * z;  
    for(int y=0; y < Ly(); ++y){    
      double ky = 2. * M_PI / Ly() * y;
      for(int x=0; x < Lx(); ++x){    
	double kx = 2. * M_PI / Lx() * x;
	std::size_t xyz_idx = xyz_to_index(x,y,z) * n_coefs();
	cx_double Rup[n_coefs()];
	calc_mat_elems(ts, delta, kx, ky, kz, Rup);
	memcpy(R_up_ + xyz_idx, Rup, sizeof(cx_double) * n_coefs());
      }
    }
  }
}

void MatElemK::get_elem(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *R) const {
  if ( is_table_set() ) {
    std::size_t xyz_idx = k_to_index(kx,ky,kz) * n_coefs();
    memcpy(R, R_up_ + xyz_idx,  sizeof(cx_double) * n_coefs());
  } else {
    calc_mat_elems(ts, delta, kx, ky, kz, R);
  }
}

MatElemK::~MatElemK(){
  if ( R_up_ != nullptr ) {
    delete[] R_up_;
  }
}


/* Member functions of MatElemN */
MatElemN::MatElemN(int Lx, int Ly, int Lz, int mu, int nu, std::vector<BondDelta> const& bonds):MatElem(Lx, Ly, Lz, 4),mu_(mu),nu_(nu),bonds_(bonds){}  // A factor of 4 comes from the combinations of 2 sublattices and 2 spins.

void MatElemN::set_table(hoppings2 const& ts, double delta){
  if ( N_ != nullptr ) {
    delete[] N_;
  }
  N_ = new cx_double[table_size()];
  build_table(ts, delta);
  assign_is_table_set(true);
}

std::size_t MatElemN::table_size() const { return Lx()*Ly()*Lz()*n_coefs(); }

void MatElemN::calc_mat_elems(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *N) const {
  /* The photon momenta are approximately set to zero. */
  calc_coef_eff_Raman_nonresonant(ts, kx, ky, kz, bonds_, mu_, nu_, N);
}

void MatElemN::build_table(hoppings2 const& ts, double delta){
  for(int z=0; z < Lz(); ++z){    
    double kz = 2. * M_PI / Lz() * z;  
    for(int y=0; y < Ly(); ++y){    
      double ky = 2. * M_PI / Ly() * y;
      for(int x=0; x < Lx(); ++x){    
	double kx = 2. * M_PI / Lx() * x;
	std::size_t xyz_idx = xyz_to_index(x,y,z) * n_coefs();
	cx_double N[n_coefs()];
	calc_mat_elems(ts, delta, kx, ky, kz, N);
	memcpy(N_ + xyz_idx, N, sizeof(cx_double) * n_coefs());
      }
    }
  }
}

void MatElemN::get_elem(hoppings2 const& ts, double delta, double kx, double ky, double kz, cx_double *N) const {
  if ( is_table_set() ) {
    std::size_t xyz_idx = k_to_index(kx,ky,kz) * n_coefs();
    memcpy(N, N_ + xyz_idx, sizeof(cx_double) * n_coefs());
  } else {
    calc_mat_elems(ts, delta, kx, ky, kz, N);
  }
}

MatElemN::~MatElemN(){
  if ( N_ != nullptr ) {
    delete[] N_;
  }
}

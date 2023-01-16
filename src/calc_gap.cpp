/*****************************************************************************
*
* Functions for calculating the gap
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "calc_gap.h"
#include "rpa_util.h"

/* Reference: A. Singh and Z. Tesanovic, PRB 41, 11457 (1990) */

bool diagonal_H(double ek1, double ek2) {  
  return ek1 * ek1 + ek2 * ek2 < 1e-12;
}
double diagonal(double ek3, double delta, double E){
  return E - ek3 - delta;
}
cx_double off_diagonal(double ek1, double ek2){
  cx_double comp(ek1, - ek2);
  return comp;  
}
cx_double calc_bk_up_in_minus(double ek1, double ek2, double ek3, double delta){
  double E = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double diag = diagonal(ek3, delta, E);
  cx_double off_diag = off_diagonal(ek1, ek2);
  double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
  return off_diag / denom;
}
cx_double calc_ak_up_in_minus(double ek1, double ek2, double ek3, double delta){
  double E = eigenenergy_HF_minus(ek1, ek2, ek3, delta);
  double diag = diagonal(ek3, delta, E);
  cx_double off_diag = off_diagonal(ek1, ek2);
  double denom = sqrt( std::norm(diag) + std::norm(off_diag) );
  return diag / denom;  
}
cx_double calc_ak_down_in_minus(double ek1, double ek2, double ek3, double delta){
  return calc_bk_up_in_minus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_down_in_minus(double ek1, double ek2, double ek3, double delta){
  return calc_ak_up_in_minus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_up_in_plus(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_bk_down_in_minus( - ek1, - ek2, ek3, delta ));  /* k + pi */
}
cx_double calc_ak_up_in_plus(double ek1, double ek2, double ek3, double delta){
  return std::conj(calc_ak_down_in_minus( - ek1, - ek2, ek3, delta ));  /* k + pi */  
}
cx_double calc_ak_down_in_plus(double ek1, double ek2, double ek3, double delta){
  return calc_bk_up_in_plus(ek1, ek2, ek3, delta);
}
cx_double calc_bk_down_in_plus(double ek1, double ek2, double ek3, double delta){
  return calc_ak_up_in_plus(ek1, ek2, ek3, delta);
}

double calc_bk_up_in(double e_free, double delta){
  double Ek = eigenenergy_HF_minus( e_free, delta );
  return  e_free / sqrt( ( Ek - delta ) * ( Ek - delta ) + e_free * e_free );
}

double calc_bk_up_out(double e_free, double delta){
  double Ek = eigenenergy_HF_plus( e_free, delta );
  
  /* Special treatment on the Fermi surface */
  double e_eps = 1e-12;
  if ( std::abs( e_free ) <= e_eps ) { return 1.; }
  else {
    return  e_free / sqrt( ( Ek - delta ) * ( Ek - delta ) + e_free * e_free );
  }
}

double calc_ak_up_in(double e_free, double delta){
  double bk = calc_bk_up_in( e_free, delta );
  return sqrt( 1. - bk * bk );
}

double calc_ak_up_out(double e_free, double delta){
  double bk = calc_bk_up_out( e_free, delta );
  return sqrt( 1. - bk * bk );
}

double calc_bk_down_in(double e_free, double delta){
  return calc_ak_up_in( e_free, delta );
}

double calc_bk_down_out(double e_free, double delta){
  return calc_ak_up_out( e_free, delta );
}

double calc_ak_down_in(double e_free, double delta){
  return calc_bk_up_in( e_free, delta );
}

double calc_ak_down_out(double e_free, double delta){
  return calc_bk_up_out( e_free, delta );
}

cx_double larger_eigenvalue(cx_double A, cx_double B, cx_double D){
  return 0.5 * ( A + D + sqrt( std::conj( A - D ) * ( A - D ) + 4. * std::conj(B) * B ) );
}

void add_to_sus_mat(cx_double& A, cx_double& B, cx_double& D, double e_free, double e_free2, double delta, cx_double omega){
  double e_eps = 1e-12;
  double E1 = eigenenergy_HF_plus( e_free2, delta );
  double E2 = eigenenergy_HF_minus( e_free, delta );
  cx_double diff_E1 = E1 - E2 + omega;
  cx_double diff_E2 = E1 - E2 - omega;
  
  double ak_up_in = calc_ak_up_in( e_free, delta );
  double ak_q_up_out = calc_ak_up_out( e_free2, delta );
  double ak_down_in = calc_ak_down_in( e_free, delta );
  double ak_q_down_out = calc_ak_down_out( e_free2, delta );
  double bk_up_in = calc_bk_up_in( e_free, delta );
  double bk_q_up_out = calc_bk_up_out( e_free2, delta );
  double bk_down_in = calc_bk_down_in( e_free, delta );
  double bk_q_down_out = calc_bk_down_out( e_free2, delta );
  
  double factor = 1.;
  
  /* Taking into account a half of the contribution from the Fermi surface */
  if ( std::abs( e_free ) <= e_eps ) {
    factor = 0.5;
  }
   
  A += factor * ( ak_up_in * ak_up_in * ak_q_down_out * ak_q_down_out / diff_E1 + ak_down_in * ak_down_in * ak_q_up_out * ak_q_up_out / diff_E2 );
  B += factor * ( ak_up_in * bk_up_in * ak_q_down_out * bk_q_down_out / diff_E1 + ak_down_in * bk_down_in * ak_q_up_out * bk_q_up_out / diff_E2 );
  D += factor * ( bk_up_in * bk_up_in * bk_q_down_out * bk_q_down_out / diff_E1 + bk_down_in * bk_down_in * bk_q_up_out * bk_q_up_out / diff_E2 );
}

void add_to_sus_mat2(hoppings const& ts, double mu, cx_double& A, cx_double& B, cx_double& C, cx_double& D, double qx, double qy, double qz, double kx, double ky, double kz, double delta, cx_double omega, bool zz){
  double eps = 1e-12;
  
  double ek1 = ts.ek1(kx, ky, kz);
  double ek2 = ts.ek2(kx, ky, kz);
  double ek3 = ts.ek3(kx, ky, kz);
  double Ek = eigenenergy_HF_minus(ek1, ek2, ek3, delta);

  cx_double ak_up = calc_ak_up_in_minus(ek1, ek2, ek3, delta);
  cx_double ak_down = calc_ak_down_in_minus(ek1, ek2, ek3, delta);
  cx_double bk_up = calc_bk_up_in_minus(ek1, ek2, ek3, delta);
  cx_double bk_down = calc_bk_down_in_minus(ek1, ek2, ek3, delta);
  
  /* Checking if the eigenenergy is below the chemical potential. */
  double mu_free = 0;  /* Assume at half filling */
  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc */
  if ( e_free > mu_free + eps ) return;
    
  /* Prefactor */
  double factor = 1.;
  
  /* On the zone boundary */
  if ( std::abs(e_free - mu_free) < eps ) {
    factor *= 0.5;
  }
  
  double diff_xm = wave_vector_in_BZ( kx - qx );
  double diff_ym = wave_vector_in_BZ( ky - qy );
  double diff_zm = wave_vector_in_BZ( kz - qz );
  
  double ek_q1m = ts.ek1(diff_xm, diff_ym, diff_zm);
  double ek_q2m = ts.ek2(diff_xm, diff_ym, diff_zm);
  double ek_q3m = ts.ek3(diff_xm, diff_ym, diff_zm);
  double Ek_qm = eigenenergy_HF_plus(ek_q1m, ek_q2m, ek_q3m, delta);
  
  cx_double diff_E2m = Ek_qm - Ek + omega;

  double diff_xp = wave_vector_in_BZ( kx + qx );
  double diff_yp = wave_vector_in_BZ( ky + qy );
  double diff_zp = wave_vector_in_BZ( kz + qz );
  
  double ek_q1p = ts.ek1(diff_xp, diff_yp, diff_zp);
  double ek_q2p = ts.ek2(diff_xp, diff_yp, diff_zp);
  double ek_q3p = ts.ek3(diff_xp, diff_yp, diff_zp);
  double Ek_qp = eigenenergy_HF_plus(ek_q1p, ek_q2p, ek_q3p, delta);
  
  cx_double diff_E1p = Ek_qp - Ek - omega;  

  /* For k-q inside the magnetic BZ, the negative sign for b comes from k -> k + pi. */
  /* For k-q outside the magnetic BZ, the negative sign for b comes from the Fourier transform. */
  double sign_A = 1.;
  double sign_B = - 1.;
  
  cx_double ak_qm_up = sign_A * calc_ak_up_in_plus( ek_q1m, ek_q2m, ek_q3m, delta);
  cx_double ak_qm_down = sign_A * calc_ak_down_in_plus( ek_q1m, ek_q2m, ek_q3m, delta);  
  cx_double bk_qm_up = sign_B * calc_bk_up_in_plus( ek_q1m, ek_q2m, ek_q3m, delta);
  cx_double bk_qm_down = sign_B * calc_bk_down_in_plus( ek_q1m, ek_q2m, ek_q3m, delta);

  cx_double ak_qp_up = sign_A * calc_ak_up_in_plus( ek_q1p, ek_q2p, ek_q3p, delta);
  cx_double ak_qp_down = sign_A * calc_ak_down_in_plus( ek_q1p, ek_q2p, ek_q3p, delta);  
  cx_double bk_qp_up = sign_B * calc_bk_up_in_plus( ek_q1p, ek_q2p, ek_q3p, delta);
  cx_double bk_qp_down = sign_B * calc_bk_down_in_plus( ek_q1p, ek_q2p, ek_q3p, delta);  
  
  if ( zz ) {
    A += factor * ( std::norm(ak_up) * std::norm(ak_qp_up) / diff_E1p + std::norm(ak_up) * std::norm(ak_qm_up) / diff_E2m );
    A += factor * ( std::norm(ak_down) * std::norm(ak_qp_down) / diff_E1p + std::norm(ak_down) * std::norm(ak_qm_down) / diff_E2m );    
    
    B += factor * ( std::conj(ak_up) * bk_up * ak_qp_up * std::conj(bk_qp_up) / diff_E1p + std::conj(bk_up) * ak_up * bk_qm_up * std::conj(ak_qm_up) / diff_E2m );
    B += factor * ( std::conj(ak_down) * bk_down * ak_qp_down * std::conj(bk_qp_down) / diff_E1p + std::conj(bk_down) * ak_down * bk_qm_down * std::conj(ak_qm_down) / diff_E2m );    
    
    C += factor * ( ak_up * std::conj(bk_up) * std::conj(ak_qp_up) * bk_qp_up / diff_E1p + bk_up * std::conj(ak_up) * std::conj(bk_qm_up) * ak_qm_up / diff_E2m );
    C += factor * ( ak_down * std::conj(bk_down) * std::conj(ak_qp_down) * bk_qp_down / diff_E1p + bk_down * std::conj(ak_down) * std::conj(bk_qm_down) * ak_qm_down / diff_E2m );    

    D += factor * ( std::norm(bk_up) * std::norm(bk_qp_up) / diff_E1p + std::norm(bk_up) * std::norm(bk_qm_up) / diff_E2m );
    D += factor * ( std::norm(bk_down) * std::norm(bk_qp_down) / diff_E1p + std::norm(bk_down) * std::norm(bk_qm_down) / diff_E2m );    
  } else {
    A += factor * ( std::norm(ak_down) * std::norm(ak_qp_up) / diff_E1p + std::norm(ak_up) * std::norm(ak_qm_down) / diff_E2m );
    
    B += factor * ( std::conj(ak_down) * bk_down * ak_qp_up * std::conj(bk_qp_up) / diff_E1p + std::conj(bk_up) * ak_up * bk_qm_down * std::conj(ak_qm_down) / diff_E2m );
    
    C += factor * ( ak_down * std::conj(bk_down) * std::conj(ak_qp_up) * bk_qp_up / diff_E1p + bk_up * std::conj(ak_up) * std::conj(bk_qm_down) * ak_qm_down / diff_E2m );

    D += factor * ( std::norm(bk_down) * std::norm(bk_qp_up) / diff_E1p + std::norm(bk_up) * std::norm(bk_qm_down) / diff_E2m );    
  }
}

/* Member functions of Polarization */
/* Polarization is depreciated: Use MatElemF instead defined in mat_elem.h. */
Polarization::Polarization(int Lx, int Ly, int Lz, int nsub):opek(){
  is_table_set_ = false;
  Lx_ = Lx;
  Ly_ = Ly;
  Lz_ = Lz;
  nsub_ = nsub;
}

void Polarization::set_q(double qx, double qy, double qz){
  qx_ = qx;
  qy_ = qy;
  qz_ = qz;
}

void Polarization::set_table(hoppings2 const& ts, double delta){
  Ppm_ = new cx_double[table_size()];
  Pzz_ = new cx_double[table_size()];  
  build_table(ts, delta);
  is_table_set_ = true;
}

long unsigned int Polarization::table_size() const { return Lx()*Ly()*Lz()*nbands()*nbands()*nsub()*nsub(); }
  
/* Assume that k = 2pi / L * m, where m is an integer. */
int Polarization::pullback(double k, int L) const {
  int l = rint(k * L / (2.*M_PI));
  while ( l < 0 ) l += L;
  while ( l >= L ) l -= L;
  /* Return 0 <= l < L */
  return l;
}
  
int Polarization::xyz_to_index(int x, int y, int z) const {
    return z * Lx() * Ly() + y * Lx() + x;
  }
  
int Polarization::k_to_index(double kx, double ky, double kz) const {    
    int lx = pullback(kx, Lx());
    int ly = pullback(ky, Ly());
    int lz = pullback(kz, Lz());    
    return xyz_to_index(lx, ly, lz);
  }

void Polarization::calc_polarization(hoppings2 const& ts, double delta, double kx, double ky, double kz, int sg1, int sg2, cx_double* ppm, cx_double* pzz) const {
  /* Adding the results */
  cx_double ek1 = ts.ek1(kx, ky, kz);
  double kx2 = kx + qx();
  double ky2 = ky + qy();
  double kz2 = kz + qz();
  cx_double ek_q1 = ts.ek1( kx2, ky2, kz2 );

  cx_double tz = ts.tz;  
  cx_vec X1up = gs_HF1(up_spin, sg1, ek1, tz, kz, delta);
  // cx_vec X1down = gs_HF1(down_spin, sg1, ek1, tz, kz, delta);
  cx_vec X2up = gs_HF1(up_spin, sg2, ek_q1, tz, kz2, delta);
  cx_vec X2down = gs_HF1(down_spin, sg2, ek_q1, tz, kz2, delta);
  
  cx_double F_zu_g1 = arma::cdot(X1up, opek.Gamma_1z * X2up);
  cx_double F_zu_gm1 = arma::cdot(X1up, opek.Gamma_2z * X2up);
	    
  /* Assume the number of sublattices is 2. */
  pzz[ 0 ] = std::norm(F_zu_g1);      
  pzz[ 1 ] = F_zu_g1 * std::conj(F_zu_gm1);
  pzz[ 2 ] = F_zu_gm1 * std::conj(F_zu_g1);
  pzz[ 3 ] = std::norm(F_zu_gm1);      

  /* <sigma_minus sigma_plus> */  
  // cx_double F_m_g1 = arma::cdot(X1down, opek.Gamma_1minus * X2up);  
  // cx_double F_m_gm1 = arma::cdot(X1down, opek.Gamma_2minus * X2up);
  // ppm[ 0 ] = std::norm(F_m_g1);  
  // ppm[ 1 ] = F_m_g1 * std::conj(F_m_gm1);
  // ppm[ 2 ] = std::conj(ppm[1]);
  // ppm[ 3 ] = std::norm(F_m_gm1);
  
  cx_double F_p_g1 = arma::cdot(X1up, opek.Gamma_1plus * X2down);
  cx_double F_p_gm1 = arma::cdot(X1up, opek.Gamma_2plus * X2down);    
  ppm[ 0 ] = std::norm(F_p_g1);      
  ppm[ 1 ] = F_p_g1 * std::conj(F_p_gm1);
  ppm[ 2 ] = F_p_gm1 * std::conj(F_p_g1);  
  ppm[ 3 ] = std::norm(F_p_gm1);
}

void Polarization::build_table(hoppings2 const& ts, double delta){
  cx_double ppm[nsub()*nsub()];
  cx_double pzz[nsub()*nsub()];
  for(int x=0; x < Lx(); x++){    
    double kx = 2. * M_PI / Lx() * x;
    for(int y=0; y < Ly(); y++){    
      double ky = 2. * M_PI / Ly() * y;
      for(int z=0; z < Lz(); z++){    
	double kz = 2. * M_PI / Lz() * z;
	long int xyz_idx = xyz_to_index(x,y,z) * nbands() * nbands() * nsub() * nsub();
	
	for(int sg1=-1; sg1<=1; sg1+=2){
	  for(int sg2=-1; sg2<=1; sg2+=2){
	    // int sg2 = - sg1; /* Opposite sign */
	    calc_polarization(ts, delta, kx, ky, kz, sg1, sg2, ppm, pzz);
	  
	    int sg1i = (sg1+1) >> 1;
	    int sg2i = (sg2+1) >> 1;
	    long int bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
	    long int xyz_bands_idx = xyz_idx + bands_idx;	  	
	    memcpy(Ppm_ + xyz_bands_idx, ppm, sizeof(cx_double)*nsub()*nsub() );
	    memcpy(Pzz_ + xyz_bands_idx, pzz, sizeof(cx_double)*nsub()*nsub() );
	  }
	}
      }
    }
  }
}

void Polarization::get_Ppm(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Ppmk) const {
  long int xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  long int bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  long int xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( Ppmk, Ppm_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}
  
void Polarization::get_Pzz(double kx, double ky, double kz, int sg1i, int sg2i, cx_double* Pzzk) const {    
  long int xyz_idx = k_to_index(kx,ky,kz) * nbands() * nbands() * nsub() * nsub();
  long int bands_idx = ((sg2i << 1) | sg1i) * nsub() * nsub();
  long int xyz_bands_idx = xyz_idx + bands_idx;    
  memcpy( Pzzk, Pzz_ + xyz_bands_idx,  sizeof(cx_double) * nsub() * nsub() );
}
  
Polarization::~Polarization(){
  if ( Ppm_ != nullptr ) {
    delete[] Ppm_;
  }

  if ( Pzz_ != nullptr ) {
    delete[] Pzz_;
  }
}

cx_double calc_prefactor_bare_res_func_bilayer(int sg1, int sg2, hoppings2 const& ts, double T, double kx, double ky, double kz, double qx, double qy, double qz, cx_double omega, double delta, double mu){
  /* k */
  cx_double ek1 = ts.ek1(kx, ky, kz);
  cx_double ek23 = ts.ek23(kx, ky, kz);
  cx_double ekz = ts.ekz(kx, ky, kz);  
  double Ek = eigenenergy_HF(sg1, ek1, ek23, ekz, ts.tz, kz, delta);
  double nk = fermi_density(Ek, kB*T, mu);
  
  /* k + q */  
  double kx2 = kx + qx;
  double ky2 = ky + qy;
  double kz2 = kz + qz;
  cx_double ek_q1 = ts.ek1( kx2, ky2, kz2 );
  cx_double ek_q23 = ts.ek23( kx2, ky2, kz2 );
  cx_double ek_qz = ts.ekz( kx2, ky2, kz2 );  
  double Ek_q = eigenenergy_HF(sg2, ek_q1, ek_q23, ek_qz, ts.tz, kz2, delta);
  double nk_q = fermi_density(Ek_q, kB*T, mu);
  
  /* Denominator */
  cx_double diff_E = omega - (Ek_q - Ek);
  
  /* Electron density at zero temperature */  
  // int n1 = ( ( sg1 + 1 ) >> 1 ) ^ 1;  /* -1 -> 1, 1 -> 0 */
  // int n2 = n1 ^ 1;
  // double n_diff = (double)(n2 - n1);
  double n_diff = nk_q - nk;
  
  cx_double prefactor = n_diff / diff_E;
  return prefactor;
}

void add_to_sus_mat4(hoppings2 const& ts, double T, double mu, arma::cx_mat& chi_pm, arma::cx_mat& chi_zz_up, arma::cx_mat& chi_zz_down, double kx, double ky, double kz, MatElemF const& me_F, double delta, cx_double omega){
// void add_to_sus_mat4(hoppings2 const& ts, double T, double mu, arma::cx_mat& chi_pm, arma::cx_mat& chi_zz_up, double kx, double ky, double kz, MatElemF const& me_F, double delta, cx_double omega){  
  
  /* Checking if the wavevector is inside the BZ. */
  double mu_free = 0;  /* Assume at half filling */
  double e_free = energy_free_electron( 1., mu_free, kx, ky );  /* ad-hoc: t=1 */
  double factor = 0.0;
  double eps = 1e-12;  
  if ( e_free > mu_free + eps ) { return; }
  else if ( std::abs(e_free - mu_free) < eps ) { factor = 0.5; } /* On the zone boundary */
  else { factor = 1.; }
      
  for(int sg1=-1; sg1<=1; sg1+=2){
    for(int sg2=-1; sg2<=1; sg2+=2){
      // int sg2 = - sg1; /* Opposite sign */
      cx_double prefactor = calc_prefactor_bare_res_func_bilayer(sg1, sg2, ts, T, kx, ky, kz, me_F.qx(), me_F.qy(), me_F.qz(), omega, delta, mu);
      prefactor *= factor;
      int sg1i = (sg1+1) >> 1;    
      int sg2i = (sg2+1) >> 1;
      cx_double F00_up[NSUBL*NSUBL];
      cx_double F00_down[NSUBL*NSUBL];            
      cx_double Fpm[NSUBL*NSUBL];
      cx_double Fzz_up[NSUBL*NSUBL];
      cx_double Fzz_down[NSUBL*NSUBL];          

      if ( me_F.is_table_set() ) {
	me_F.get_pm(kx, ky, kz, sg1i, sg2i, Fpm);
	me_F.get_zz_up(kx, ky, kz, sg1i, sg2i, Fzz_up);
	me_F.get_zz_down(kx, ky, kz, sg1i, sg2i, Fzz_down);	
      } else {
	me_F.calc_mat_elems(ts, delta, kx, ky, kz, sg1, sg2, F00_up, F00_down, Fpm, Fzz_up, Fzz_down);
      }

      chi_pm(0,0) += prefactor * Fpm[0];
      chi_pm(0,1) += prefactor * Fpm[1];
      chi_pm(1,0) += prefactor * Fpm[2];
      chi_pm(1,1) += prefactor * Fpm[3];
    
      chi_zz_up(0,0) += prefactor * Fzz_up[0];
      chi_zz_up(0,1) += prefactor * Fzz_up[1];
      chi_zz_up(1,0) += prefactor * Fzz_up[2];
      chi_zz_up(1,1) += prefactor * Fzz_up[3];

      chi_zz_down(0,0) += prefactor * Fzz_down[0];
      chi_zz_down(0,1) += prefactor * Fzz_down[1];
      chi_zz_down(1,0) += prefactor * Fzz_down[2];
      chi_zz_down(1,1) += prefactor * Fzz_down[3];      
    }
  }
}

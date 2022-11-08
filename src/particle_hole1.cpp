/*****************************************************************************
*
* Functions for calculating the one particle-hole-pair states.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include "array3.h"
#include "particle_hole1.h"
#include "calc_gap.h"

extern int valence_index;
extern int conduction_index;

typedef rpa::array3 vec3;

std::tuple<std::size_t, int> index_to_kidx_spin(std::size_t idx){
  int spin = idx & 1;
  std::size_t kidx = (idx >> 1);
  return std::make_tuple(kidx, spin);
}

std::size_t kidx_spin_to_index(std::size_t ki, int si){
  return (ki << 1) | si;  
}

int sigma_z(int sigma1, int sigma2){
  if ( sigma1 == sigma2 ) {
    if ( sigma1 == up_spin ) { return 1; }
    else { return - 1; }
  } else {
    return 0;
  }
}
      
double calc_factor_V(double U, int s1, int s2, int s3, int s4){
  int factor = 1;
      
  /* sigma^z_{s1 s4} */
  factor *= sigma_z(s1, s4);

  /* sigma^z_{s2 s3} */
  factor *= sigma_z(s2, s3);      
      
  return 0.5 * U * factor;
  // return - 0.5 * U * factor;
  
}

void pop_back_operator(std::vector<ElecIndex>& eis, std::vector<std::vector<cx_double>>& mat_elems){
  eis.pop_back();
  for(int g=0; g < NSUBL; ++g){
    mat_elems[g].pop_back();
  }
}

int k_to_index(double k, std::size_t L){
  if ( L == 1 ) { return 0; }
  else {
    double k2 = wave_vector_in_BZ(k);
    double k1 = 2 * M_PI / L;
    return (int)(std::round((k2 + M_PI) / k1));
  }
}

std::size_t k_to_index(double kx, double ky, double kz, std::size_t Lx, std::size_t Ly, std::size_t Lz){
  int xi = k_to_index(kx, Lx);
  int yi = k_to_index(ky, Ly);
  int zi = k_to_index(kz, Lz);
  return zi * Lx * Ly + yi * Lx + xi;  
}

std::size_t k_to_index(vec3 kvec, std::size_t Lx, std::size_t Ly, std::size_t Lz){
  return k_to_index(kvec[0], kvec[1], kvec[2], Lx, Ly, Lz);
}

std::vector<double> make_1d_k_vector(int L){
  if ( L == 1 ) {
    std::vector<double> ks;
    ks.push_back(0.);
    return ks;
  } else {
    double k1 = 2.0 * M_PI / L;
    std::vector<double> ks;
    for(int x = - L/2; x < L/2; ++x){
      ks.push_back(k1 * x);
    }
    return ks;
  }
}

void make_wave_vector_tables(rpa::parameters const& pr, std::vector<std::size_t>& k_index_table, std::vector<vec3>& ks){
  /* Parameters */
  int Lx = pr.Lx;
  int Ly = pr.Ly;
  int Lz = pr.Lz;
  int N = Lx * Ly * Lz;
  
  /* Assume the ordering wave vector to be (pi,pi,pi). */
  double Qx = Lx == 1 ? 0. : M_PI;
  double Qy = Ly == 1 ? 0. : M_PI;
  double Qz = Lz == 1 ? 0. : M_PI;  
  vec3 Qvec({Qx, Qy, Qz});

  /* Wave vectors for each axis */
  std::vector<double> kxs = make_1d_k_vector(Lx);
  std::vector<double> kys = make_1d_k_vector(Ly);
  std::vector<double> kzs = make_1d_k_vector(Lz);
  
  k_index_table.clear();
  k_index_table.resize(N);
  
  ks.clear();
  
  std::vector<bool> is_assigned(N, false);
  
  /* Building ks. */
  for(double kz: kzs){
    for(double ky: kys){
      for(double kx: kxs){	
	/* Checking if the wavevector is inside the BZ. */
	double factor = BZ_factor(kx, ky);
	
	if ( std::abs(factor) < 1e-12 ) { continue; }
	
	// if ( std::abs(factor - 0.5) < 1e-12 && ky > - 1e-12 ) { continue; }
	if ( std::abs(factor - 0.5) < 1e-12 && kx > - 1e-12 ) { continue; }

	std::size_t k_idx = k_to_index(kx, ky, kz, Lx, Ly, Lz);
	k_index_table[k_idx] = ks.size();
	is_assigned[k_idx] = true;
	ks.push_back(vec3({kx, ky, kz}));
      }
    }
  }

  /* Completing the index table for the unfolded reciprocal space. */
  for(double kz: kzs){
    for(double ky: kys){
      for(double kx: kxs){  
	std::size_t k_idx = k_to_index(kx, ky, kz, Lx, Ly, Lz);
	if ( !is_assigned[k_idx] ) {
	  /* Adding Qvec */
	  std::size_t k_idx2 = k_to_index(kx + Qvec[0], ky + Qvec[1], kz + Qvec[2], Lx, Ly, Lz);
	  assert(is_assigned[k_idx2]);
	  k_index_table[k_idx] = k_index_table[k_idx2];
	  is_assigned[k_idx] = true;	  	  
	}
      }
    }
  }
}

void calc_particle_hole1(rpa::parameters const& pr_, hoppings2 const& ts, double delta, vec& vals, cx_mat& Uph1){
  /* Parameters */
  rpa::parameters pr(pr_);  
  int Lx = pr.Lx;
  int Ly = pr.Ly;
  int Lz = pr.Lz;
  int N = Lx * Ly * Lz;
  double U = pr.U;
    
  /* Index of wave vectors*/  
  std::vector<std::size_t> k_index_table;
  std::vector<vec3> ks;
  make_wave_vector_tables(pr, k_index_table, ks);
  std::size_t nk = ks.size();
  assert(nk == N / NSUBL);  
  
  /* Unitary matrix */
  std::vector<cx_mat> Uks;
  for(std::size_t ki=0; ki < nk; ++ki){
    /* Wave vector */
    double kx = ks[ki][0];
    double ky = ks[ki][1];
    double kz = ks[ki][2];
    cx_double ek1 = ts.ek1(kx, ky, kz);
    cx_mat Uk = gs_HF(ek1, ts.tz, kz, delta);
    Uks.push_back(Uk);
  }
  
  std::size_t dim = nk * 2;    
  cx_mat H1(dim, dim, arma::fill::zeros);
  
  /* Building the Hamiltonian. */
  for(std::size_t i=0; i < dim; ++i){
    if ( true ) {          
      std::size_t ki = 0;
      int sigmai_idx = 0, sigmai = 0, sigmai_bar = 0;
      std::tie(ki, sigmai_idx) = index_to_kidx_spin(i);
      sigmai = index_to_spin(sigmai_idx);
      sigmai_bar = sigmai == up_spin ? down_spin : up_spin;
      cx_mat Ui = Uks[ki];
      double kx = ks[ki][0];
      double ky = ks[ki][1];
      double kz = ks[ki][2];
    
      /* Eigenenergy */
      cx_double ek1 = ts.ek1(kx, ky, kz);
      cx_double ek23 = ts.ek23(kx, ky, kz);
      cx_double ekz = ts.ekz(kx, ky, kz);
      double ek_plus = eigenenergy_HF(1., ek1, ek23, ekz, ts.tz, kz, delta);
      double ek_minus = eigenenergy_HF(-1., ek1, ek23, ekz, ts.tz, kz, delta);
      double e_diff = ek_plus - ek_minus;
    
      /* Diagonal term */
      H1(i,i) += e_diff;
    } /* end if */
    
    /* Band index */
    std::vector<int> valence({valence_index});
    std::vector<int> conduction({conduction_index});
    
    /* Electronic state */
    ElecState esi(valence, conduction);
    
    /* Taking the sigma z form into account. */
    auto valid_spins = [](int sigma1, int sigma2, int sigma3, int sigma4){
      // if ( sigma1 == sigma4 && sigma2 == sigma3 && sigma1 != sigma2 ) { return true; }
      // else { return false; }

      if ( sigma1 == sigma4 && sigma2 == sigma3 ) { return true; }
      else { return false; }      
    };

    auto valid_combination = [](std::size_t kki, std::size_t kpi, std::pair<int, int> bs1, std::pair<int, int> bs2, std::pair<int, int> bs3, std::pair<int, int> bs4){
      // /* (v c v c), (c v c v), (v c c v), or (c v v c) */	
      // if ( bs1.first != bs2.first && bs3.first != bs4.first ) {}
      // else { return false; }
      
      // /* (up down down up) or (down up up down) */      
      // if ( bs1.second == bs4.second && bs2.second == bs3.second && bs1.second != bs2.second ) {}
      // else { return false; }

      /* sigma1 == sigma4 and sigma2 == sigma3 */
      if (bs1.second == bs4.second && bs2.second == bs3.second) {}
      else { return false; }          

      return true;
    };    

    /* Setting single-particle-hole pair state. */
    std::size_t ki = 0;
    int sigmai_idx = 0;
    std::tie(ki, sigmai_idx) = index_to_kidx_spin(i);
    bool test = true;    
    test = esi.elec_creation(ElecIndex(ki, conduction_index, sigmai_idx));
    test = esi.elec_annihilation(ElecIndex(ki, valence_index, sigmai_idx));
    
    std::vector<std::pair<int, int>> band_spins;
    band_spins.push_back(std::pair(0, 0));
    band_spins.push_back(std::pair(0, 1));
    band_spins.push_back(std::pair(1, 0));
    band_spins.push_back(std::pair(1, 1));
    std::vector<bool> is_elec_creation = {true, true, false, false};
    for(std::size_t kki = 0; kki < nk; ++kki){
      vec3 kk = ks[kki];
      for(std::size_t kpi = 0; kpi < nk; ++kpi){
	vec3 kp = ks[kpi];

	/* kq = 0 */
	// vec3 kq({0, 0, 0});	
	// vec3 kq({M_PI, M_PI, M_PI});
	
	for(std::size_t kqi = 0; kqi < nk; ++kqi){
	  vec3 kq = ks[kqi];		  

	  vec3 kk_plus_q = kk + kq;
	  std::size_t kk_plus_qi = k_index_table[k_to_index(kk_plus_q, Lx, Ly, Lz)];

	  vec3 kp_minus_q = kp - kq;
	  std::size_t kp_minus_qi = k_index_table[k_to_index(kp_minus_q, Lx, Ly, Lz)];

	  std::unordered_set<std::size_t> kset = {kki, kpi, kk_plus_qi, kp_minus_qi};

	  /* Number of distinct wave vectors <= 2 */
	  if ( kset.size() <= 2 ) {	      
	    std::vector<ElecIndex> eis;
	    std::vector<std::vector<cx_double>> mat_elems(NSUBL);  // for each gamma and operator index
		    
	    for(std::pair<int, int> bs1: band_spins){
	      /* Operator 1 */
	      std::size_t ok1 = kki;
	      ElecIndex ei1(ok1, bs1.first, bs1.second);
	      eis.push_back(ei1);

	      /* Matrix element 1 */	    
	      int sigma1 = index_to_spin(bs1.second);
	      int column_index = band_spin_index(bs1.first, sigma1);	    
	      for(int g=0; g < NSUBL; ++g){
		int row_index = sublattice_spin_index(g, sigma1);
		cx_double Uelem = Uks[ok1](row_index, column_index);
		if ( is_elec_creation[0] ) {
		  mat_elems[g].push_back(std::conj(Uelem));
		} else {
		  mat_elems[g].push_back(Uelem);
		}
	      }  /* end for g */
		    
	      for(std::pair<int, int> bs2: band_spins){
		/* Operator 2 */	      
		std::size_t ok2 = kpi;
		ElecIndex ei2(ok2, bs2.first, bs2.second);
		eis.push_back(ei2);

		/* Matrix element 2 */	      
		int sigma2 = index_to_spin(bs2.second);
		int column_index = band_spin_index(bs2.first, sigma2);	      
		for(int g=0; g < NSUBL; ++g){
		  int row_index = sublattice_spin_index(g, sigma2);
		  cx_double Uelem = Uks[ok2](row_index, column_index);
		  if ( is_elec_creation[1] ) {
		    mat_elems[g].push_back(std::conj(Uelem));
		  } else {
		    mat_elems[g].push_back(Uelem);
		  }
		}  /* end for g */
		    
		for(std::pair<int, int> bs3: band_spins){
		  /* Operator 3 */
		  std::size_t ok3 = kp_minus_qi;
		  ElecIndex ei3(ok3, bs3.first, bs3.second);
		  eis.push_back(ei3);

		  /* Matrix element 3 */	      
		  int sigma3 = index_to_spin(bs3.second);
		  int column_index = band_spin_index(bs3.first, sigma3);		
		  for(int g=0; g < NSUBL; ++g){
		    int row_index = sublattice_spin_index(g, sigma3);
		    cx_double Uelem = Uks[ok3](row_index, column_index);
		    if ( is_elec_creation[2] ) {
		      mat_elems[g].push_back(std::conj(Uelem));
		    } else {
		      mat_elems[g].push_back(Uelem);
		    }		  
		  }  /* end for g */
		    
		  for(std::pair<int, int> bs4: band_spins){
		    /* Checking if the sequence of the spin indices is allowed. */
		    if ( valid_combination(kki, kpi, bs1, bs2, bs3, bs4) ) {
		      /* Operator 4 */
		      std::size_t ok4 = kk_plus_qi;
		      ElecIndex ei4(ok4, bs4.first, bs4.second);
		      eis.push_back(ei4);
		    
		      /* Matrix element 4 */	      
		      int sigma4 = index_to_spin(bs4.second);
		      int column_index = band_spin_index(bs4.first, sigma4);
		      for(int g=0; g < NSUBL; ++g){
			int row_index = sublattice_spin_index(g, sigma4);
			cx_double Uelem = Uks[ok4](row_index, column_index);
		    
			if ( is_elec_creation[3] ) {
			  mat_elems[g].push_back(std::conj(Uelem));
			} else {
			  mat_elems[g].push_back(Uelem);
			}		      
		      }  /* end for g */
		    
		      assert(eis.size() == is_elec_creation.size());

		      /* Applying the operators. */
		      ElecState esf(esi);
		      bool valid = true;
		    
		      /* Backward */
		      for(int opi=eis.size()-1; opi >= 0; --opi){		    
			if ( is_elec_creation[opi] ) {
			  valid = esf.elec_creation(eis[opi]);
			} else {
			  valid = esf.elec_annihilation(eis[opi]);
			}
			if ( !valid ) { break; }	
		      } /* end for opi */
		    
		      if ( valid ) {
			if ( esf.particle_hole_pair_state() ) {
			  int n_particles = esf.num_particles();
			  std::size_t j = 0;
			  if ( n_particles <= 1 ) {
			    cx_double sum = 0;
			    for(int g=0; g < NSUBL; ++g){
			      cx_double mat_elem = 1.0;
			      for(cx_double elem: mat_elems[g]){
				mat_elem *= elem;
			      }
			      sum += mat_elem;
			    }  /* end for g */
			    sum *= U / (double)N;
			  
			    if ( n_particles == 0 ) {
			      j = 0;
			    } else if ( n_particles == 1 ) {
			      auto eif_e = esf.particles.begin();
			      auto eif_h = esf.holes.begin();
			      
			      /* Zero total momentum and Sz=0 */
			      assert(eif_e->k_index == eif_h->k_index && eif_e->spin_index == eif_h->spin_index);
			    
			      j = kidx_spin_to_index(eif_e->k_index, eif_e->spin_index);			  
			    }

			    if ( n_particles == 1 ) {
			      H1(j,i) += sum;

			      // std::cout << i << "  " << j << "  " << sum << "  " << ok1 << "  " << bs1.first << " " << bs1.second << "    " << ok2 << "  " << bs2.first << " " << bs2.second << "    " << ok3 << "  " << bs3.first << " " << bs3.second << "    " << ok4 << "  " << bs4.first << " " << bs4.second << " " << std::endl;
			    }
			  }
			}
		      }
		      pop_back_operator(eis, mat_elems);
		    }  /* end if */
		  }  /* end for bs4 */
		  pop_back_operator(eis, mat_elems);
		}  /* end for bs3 */
		pop_back_operator(eis, mat_elems);
	      }  /* end for bs2 */
	      pop_back_operator(eis, mat_elems);
	    }  /* end for bs1 */
	  } /* end if */
	}  /* end for kqi */
      }  /* end for kpi */
    }  /* end for kki */
  }

  // // for check
  // std::cout << H1 << std::endl;
  
  /* Diagonalization */
  arma::eig_sym(vals, Uph1, H1);
}

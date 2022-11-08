/*****************************************************************************
*
* Functions for calculating the wave function
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <complex>
#ifdef WITH_OpenMP
#include <omp.h>
#endif
#include "calc_intensity.h"
#include "calc_spectrum.h"
#include "self_consistent_eq.h"
#include "calc_gap.h"
#include "calc_chemical_potential.h"
#include "calc_wave_func.h"
#include "BinarySearch.h"
#include "particle_hole1.h"
#include "calc_two_site.h"
#include "calc_gap_two_site.h"

extern int valence_index;
extern int conduction_index;

void calc_particle_hole_two_site(rpa::parameters const& pr_, hoppings2 const& ts, double delta, vec& vals, cx_mat& Uph1){
  rpa::parameters pr(pr_);
  pr.Lx = 2;
  pr.Ly = 1;
  pr.Lz = 1;
  int N = 2;
  
  /* Unitary matrix */
  cx_mat Uk = gs_HF(ts.ek1(0, 0, 0), 0., 0., delta);

  std::size_t dim = 4;
  cx_mat H1(dim, dim, arma::fill::zeros);

  /* Diagonal term */
  H1(0,0) += 2.0 * ts.ek1(0, 0, 0);      // Ground state
  H1(1,1) += 2.0 * ts.ek1(M_PI, 0, 0);   // Two particle-hole pairs
  H1(2,2) += 0.0;     // One particle-hole pair
  H1(3,3) += 0.0;     // One particle-hole pair  
			   
  for(int i=0; i < dim; ++i){
    /* Band index */			   
    std::vector<int> valence({valence_index});
    std::vector<int> conduction({conduction_index});
    
    /* Electronic state */
    ElecState esi(valence, conduction);
    
    /* Taking the sigma z form into account. */
    auto valid_spins = [](int sigma1, int sigma2, int sigma3, int sigma4){
      if ( sigma1 == sigma4 && sigma2 == sigma3 && sigma1 != sigma2 ) { return true; }
      else { return false; }
    };

    auto valid_combination = [](std::pair<int, int> bs1, std::pair<int, int> bs2, std::pair<int, int> bs3, std::pair<int, int> bs4){
      // if ( false ) {      
      // if ( true ) {
      // if ( kki == kpi ) {
	
      // 	/* (v c v c), (c v c v), (v c c v), or (c v v c) */	
      // 	if ( bs1.first != bs2.first && bs3.first != bs4.first ) {}
      // 	else { return false; }
      // } else {
      // 	/* (v c v c) or (c v c v) */	
      // 	if ( bs1.first == bs3.first && bs2.first == bs4.first && bs1.first != bs2.first ) {}
      // 	else { return false; }
      // }

      /* (up down down up) or (down up up down) */      
      if ( bs1.second == bs4.second && bs2.second == bs3.second ) {}
      // if ( bs1.second == bs4.second && bs2.second == bs3.second && bs1.second != bs2.second ) {}      
      else { return false; }

      return true;
    };

    
    // for check
    if ( i == 0 ) {
      /* Ground state */
    } else if ( i == 1 ) {
      /* Two particle-hole pairs */      
      bool test = true;
      test = esi.elec_creation(ElecIndex(0, conduction_index, 0));
      test = esi.elec_annihilation(ElecIndex(0, valence_index, 0));      
      test = esi.elec_creation(ElecIndex(0, conduction_index, 1));
      test = esi.elec_annihilation(ElecIndex(0, valence_index, 1));            
    } else if ( i == 2 ) {
      /* Up-spin particle-hole pair */            
      bool test = true;      
      test = esi.elec_creation(ElecIndex(0, conduction_index, 0));      
      test = esi.elec_annihilation(ElecIndex(0, valence_index, 0));
    } else if ( i == 3 ) {
      /* Down-spin particle-hole pair */                  
      bool test = true;      
      test = esi.elec_creation(ElecIndex(0, conduction_index, 1));      
      test = esi.elec_annihilation(ElecIndex(0, valence_index, 1));
    } else {}
    
    std::vector<std::pair<int, int>> band_spins;
    band_spins.push_back(std::pair(0, 0));
    band_spins.push_back(std::pair(0, 1));
    band_spins.push_back(std::pair(1, 0));
    band_spins.push_back(std::pair(1, 1));
    std::vector<bool> is_elec_creation = {true, true, false, false};
	      
    std::vector<ElecIndex> eis;
    std::vector<std::vector<cx_double>> mat_elems(NSUBL);  // for each gamma and operator index
		    
    for(std::pair<int, int> bs1: band_spins){
      /* Operator 1 */
      ElecIndex ei1(0, bs1.first, bs1.second);
      eis.push_back(ei1);

      /* Matrix element 1 */	    
      int sigma1 = index_to_spin(bs1.second);
      int column_index = band_spin_index(bs1.first, sigma1);	    
      for(int g=0; g < NSUBL; ++g){
	int row_index = sublattice_spin_index(g, sigma1);
	cx_double Uelem = Uk(row_index, column_index);
	if ( is_elec_creation[0] ) {
	  mat_elems[g].push_back(std::conj(Uelem));
	} else {
	  mat_elems[g].push_back(Uelem);
	}
      }  /* end for g */
		    
      for(std::pair<int, int> bs2: band_spins){
	/* Operator 2 */	      
	ElecIndex ei2(0, bs2.first, bs2.second);
	eis.push_back(ei2);

	/* Matrix element 2 */	      
	int sigma2 = index_to_spin(bs2.second);
	int column_index = band_spin_index(bs2.first, sigma2);	      
	for(int g=0; g < NSUBL; ++g){
	  int row_index = sublattice_spin_index(g, sigma2);
	  cx_double Uelem = Uk(row_index, column_index);
	  if ( is_elec_creation[1] ) {
	    mat_elems[g].push_back(std::conj(Uelem));
	  } else {
	    mat_elems[g].push_back(Uelem);
	  }
	}  /* end for g */
		    
	for(std::pair<int, int> bs3: band_spins){
	  /* Operator 3 */
	  ElecIndex ei3(0, bs3.first, bs3.second);
	  eis.push_back(ei3);

	  /* Matrix element 3 */	      
	  int sigma3 = index_to_spin(bs3.second);
	  int column_index = band_spin_index(bs3.first, sigma3);		
	  for(int g=0; g < NSUBL; ++g){
	    int row_index = sublattice_spin_index(g, sigma3);
	    cx_double Uelem = Uk(row_index, column_index);
	    if ( is_elec_creation[2] ) {
	      mat_elems[g].push_back(std::conj(Uelem));
	    } else {
	      mat_elems[g].push_back(Uelem);
	    }		  
	  }  /* end for g */
		    
	  for(std::pair<int, int> bs4: band_spins){
	    /* Checking if the sequence of the spin indices is allowed. */

	    // double factor_V = calc_factor_V(U, bs1.second, bs2.second, bs3.second, bs4.second);
	    // if ( std::abs(factor_V) > 1e-12 && valid_combination(kki, kpi, bs1, bs2, bs3, bs4) ) {
	    if ( valid_combination(bs1, bs2, bs3, bs4) ) {
	      // if ( valid_combination(bs1, bs2, bs3, bs4) ) {		      
	      // if ( valid_spins(bs1.second, bs2.second, bs3.second, bs4.second) ) {
		      
	      /* Operator 4 */
	      ElecIndex ei4(0., bs4.first, bs4.second);
	      eis.push_back(ei4);
		    
	      /* Matrix element 4 */	      
	      int sigma4 = index_to_spin(bs4.second);
	      int column_index = band_spin_index(bs4.first, sigma4);
	      for(int g=0; g < NSUBL; ++g){
		int row_index = sublattice_spin_index(g, sigma4);
		cx_double Uelem = Uk(row_index, column_index);
		    
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

		// // for check
		// std::cerr << i << "  valid:  n_e = " << esf.particles.size() << "  n_h = " << esf.holes.size() << std::endl;
		      
		if ( esf.particle_hole_pair_state() ) {

		  // // for check
		  // std::cerr << "particle hole pair" << std::endl;
		      
		  cx_double sum = 0;
		  for(int g=0; g < NSUBL; ++g){
		    cx_double mat_elem = 1.0;
		    for(cx_double elem: mat_elems[g]){
		      mat_elem *= elem;

		      // for check
		      std::cout << "  " << elem;
		    }
		    sum += mat_elem;

		    // for check
		    std::cout << sum << std::endl;
				
		  }  /* end for g */
		  // sum *= factor_V / (double)N;
		  sum *= pr.U / (double)N;

		  // auto eif_e = esf.particles.begin();
		  // auto eif_h = esf.holes.begin();

		  // /* Zero total momentum and Sz=0 */
		  // assert(eif_e->k_index == eif_h->k_index && eif_e->spin_index == eif_h->spin_index);

		  int n_particles = esf.num_particles();
		  int j = 0;		  
		  if ( n_particles == 0 ) {
		    j = 0;
		  } else if ( n_particles == 2 ) {
		    j = 1;			      
		  } else if ( n_particles == 1 ) {
		    auto eif_e = esf.particles.begin();
		    if ( eif_e->spin_index == 0 ) {
		      j = 2;
		    } else {
		      j = 3;
		    }
		  } else {}
			    
		  // // for check
		  // std::cerr << "one pair" << std::endl;
		      

		  // // for check
		  // std::cout << eif_e->k_index << "  " << eif_h->k_index << "  " << eif_e->spin_index << "  " << eif_h->spin_index << std::endl;
			      
		  H1(j,i) += sum;

		  std::cout << i << "  " << j << "  " << sum << "  " << "  " << bs1.first << " " << bs1.second << "    " << "  " << bs2.first << " " << bs2.second << "    " << "  " << bs3.first << " " << bs3.second << "    " << "  " << bs4.first << " " << bs4.second << " " << std::endl;		

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

		    // std::cerr << __FILE__ << "  " << __LINE__ << std::endl;

  // for check
  std::cout << H1 << std::endl;
  
  /* Diagonalization */
  arma::eig_sym(vals, Uph1, H1);

}

void calc_two_site(path& base_dir, rpa::parameters const& pr){
  /* Getting parameters */
  double U = pr.U;
  double T = pr.T;

  /* Hopping parameters */
  std::unique_ptr<hoppings_two_site> ts = hoppings_two_site::mk_two_site(pr);

  /* Calculate the chemical potential and the charge gap. */
  double mu = 0;
  double delta = solve_self_consistent_eq_two_site(*ts, U, T);  
  std::cout << "delta = " << delta << std::endl;  

  /* Output */
  ofstream out_gap;
  out_gap.open(base_dir / "gap.text");
  out_gap << "mu = " << mu << std::endl;
  
  /* Polarization */
  int Lx = 2;
  int Ly = 1;
  int Lz = 1;  
  Polarization Pz( Lx, Ly, Lz, NSUBL );

  auto calc_gaps = [&](double q){
    /* Setting the ordering vector */
    Pz.set_q(q, 0, 0);
    
    /* Calculating gaps */
    double omega_T = 0, omega_L = 0, omega_ph = 0;
    bool return_upper = true;
    // bool verbose = true;
    bool verbose = false;    
    std::tie(omega_T, omega_L, omega_ph) = calc_gap_two_site(*ts, mu, U, T, delta, Pz, return_upper, verbose);
    
    return std::make_tuple(omega_T, omega_L, omega_ph);
  };

  /* Calculating gaps */
  double omega_T = 0, omega_L = 0, omega_ph = 0;
  std::tie(omega_T, omega_L, omega_ph) = calc_gaps(0);
  // std::tie(omega_T, omega_L, omega_ph) = calc_gaps(M_PI);
  
  /* Output */
  out_gap << "omega_T = " << omega_T << std::endl;
  out_gap << "omega_L = " << omega_L << std::endl;
  out_gap << "omega_ph = " << omega_ph << std::endl;

  // /* Precision */
  // int prec = 10;
  // int pw = prec + 10;

  /* Calculating the scattering states of a single particle-hole pair. */
  vec vals;
  cx_mat Uph1;
  calc_particle_hole_two_site(pr, *ts, delta, vals, Uph1);

  // for check
  std::cout << "eigval" << std::endl;
  std::cout << vals << std::endl;
  
  /* Closing the output files. */
  out_gap.close();
}
			 

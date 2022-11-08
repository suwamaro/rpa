/*****************************************************************************
*
* Functions for calculating the one particle-hole-pair states.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#ifndef _PARTICLE_HOLE1_H
#define _PARTICLE_HOLE1_H

#include <unordered_set>
#include <boost/functional/hash.hpp>
#include "rpa.h"
#include "rpa_util.h"

class ElecIndex {
public:
  ElecIndex():k_index(0), band_index(0), spin_index(0){}
  explicit ElecIndex(std::size_t k, int b, int s):k_index(k), band_index(b), spin_index(s){}
  std::size_t k_index;
  int band_index;
  int spin_index;
  
  bool operator ==(ElecIndex const& ei) const {
    if ( this->k_index == ei.k_index
	 && this->band_index == ei.band_index
	 && this->spin_index == ei.spin_index ) { return true; }
    else { return false; }
  }
};

class ElecIndexHasher {
public:
  std::size_t operator()(ElecIndex const& ei) const {
    std::size_t seed = 0;
    boost::hash_combine(seed, ei.k_index);
    boost::hash_combine(seed, ei.band_index);
    boost::hash_combine(seed, ei.spin_index);
    return seed;
  }
};

class ElecIndexEqual {
public:
  bool operator()(ElecIndex const& lhs, ElecIndex const& rhs) const {
    return lhs == rhs;
  }
};
  
class ElecState {
public:
  ElecState(){}
  explicit ElecState(std::vector<int> const& v, std::vector<int> const& c): valence(v.begin(), v.end()), conduction(c.begin(), c.end()) {}
  explicit ElecState(ElecState const& es){
    valence = es.valence;
    conduction = es.conduction;
    particles = es.particles;
    holes = es.holes;
  }
  std::unordered_set<int> valence;
  std::unordered_set<int> conduction;  
  std::unordered_set<ElecIndex, ElecIndexHasher, ElecIndexEqual> particles;
  std::unordered_set<ElecIndex, ElecIndexHasher, ElecIndexEqual> holes;

  bool particle_hole_pair_state() const {
    return particles.size() == holes.size();
  }
  
  int num_particles() const {
    return particles.size();
  }
  
  bool elec_creation(ElecIndex const& ei){
    if (conduction.contains(ei.band_index)) {
      if (particles.contains(ei)) { return false; }
      else {
	particles.insert(ei);
	return true;
      }
    } else if (valence.contains(ei.band_index)) {
      if (holes.contains(ei)) {
	holes.erase(ei);
	return true;
      } else {
	return false;
      }
    } else {
      std::cerr << "Band index " << ei.band_index << " is not labeled as a valence or conduction band in function " << __func__ << ".\n";
      std::exit(EXIT_FAILURE);
    }
  }

  bool elec_annihilation(ElecIndex const& ei){
    if (conduction.contains(ei.band_index)) {
      if (particles.contains(ei)) {
	particles.erase(ei);
	return true;
      } else {
	return false;
      }
    } else if (valence.contains(ei.band_index)) {
      if (holes.contains(ei)){
	return false;
      } else {
	holes.insert(ei);
	return true;
      }
    } else {
      std::cerr << "Band index " << ei.band_index << " is not labeled as a valence or conduction band in function " << __func__ << ".\n";
      std::exit(EXIT_FAILURE);
    }
  }
};

void pop_back_operator(std::vector<ElecIndex>& eis, std::vector<std::vector<cx_double>>& mat_elems);
  
void calc_particle_hole1(rpa::parameters const& pr, hoppings2 const& ts, double delta, vec& vals, cx_mat& Uph1);

#endif // _PARTICLE_HOLE1_H

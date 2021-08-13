/*****************************************************************************
*
* Functions for mean field calculations.
*
* Copyright (C) 2018 by Hidemaro Suwa
* e-mail:suwamaro@phys.s.u-tokyo.ac.jp
*
*****************************************************************************/

#include <cuba.h>
#include "rpa.h"
#include "cuba_helper.h"
#include "hoppings.h"

class MeanField {
public:
  explicit MeanField(std::size_t n);
  virtual void set_eigen_energies() = 0;
  virtual void set_eigen_vectors() = 0;    
  double E(std::size_t i) const { return E_(i); }
  const cx_mat *U() const { return &U_; }
  // const cx_mat *U_dagger() const { return &U_dagger_; }  
  cx_double U(std::size_t i, std::size_t j) const { return U(i,j); }
  // cx_double U_dagger(std::size_t i, std::size_t j) const { return U_dagger(i,j); }
  virtual ~MeanField(){}

protected:
  vec E_;  
  cx_mat U_;
  // cx_mat U_dagger_;
};

class MeanFieldBilayer : public MeanField {
public:
  explicit MeanFieldBilayer();  
  void set_eigen_energies();
  void set_eigen_vectors();
  void set_parameters(hoppings_bilayer2 const& h, double delta, double* kvec);    
  cx_vec calc_eigen_vector(int sign, int spin, cx_double xk, double zk) const;

private:
  hoppings_bilayer2 hb_;
  double delta_;
  double kvec_[3];
};

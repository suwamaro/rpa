/***************************************************************************
* 3d vector (array)
*
* Copyright (C) 2018 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/
#include <cmath>
#include "array3.h"

namespace rpa {
  /* Member functions of array3 */
  array3::array3() { clear(); }
  array3::array3(std::size_t i) { clear(); }
  array3::array3(double x, double y, double z) {
    super_type::operator[](0) = x;
    super_type::operator[](1) = y;
    super_type::operator[](2) = z;
  }
  array3::array3(array3 const& v){
    for (int i = 0; i < 3; ++i) super_type::operator[](i) = v[i];
  }
  array3::array3(std::initializer_list<double> init){
    assert(init.end() - init.begin() == 3);
    std::copy(init.begin(), init.end(), super_type::begin());
  }
  void array3::clear() { for (int i = 0; i < 3; ++i) super_type::operator[](i) = 0; }
  double array3::norm2() const { return (*this) * (*this); }
  double array3::norm() const { return sqrt(norm2()); }
  double array3::operator*(array3 const& v) const {
    double res = 0;
    for (int i = 0; i < 3; ++i) res += super_type::operator[](i) * v[i];
    return res;
  }
  array3 array3::operator*=(double b) {
    for (int i = 0; i < 3; ++i) super_type::operator[](i) *= b;
    return *this;
  }
  array3 array3::operator/=(double b) {
    for (int i = 0; i < 3; ++i) super_type::operator[](i) /= b;
    return *this;
  }    
  array3& array3::operator+=(array3 const& v) {
    for (int i = 0; i < 3; ++i) super_type::operator[](i) += v[i];
    return *this;
  }
  array3& array3::operator-=(array3 const& v) {
    for (int i = 0; i < 3; ++i) super_type::operator[](i) -= v[i];
    return *this;
  }
  array3 array3::operator-() const {
    array3 v2;
    v2 -= *this;
    return v2;
  }
  array3& array3::operator=(array3 const& v) {
    for (int i = 0; i < 3; ++i) super_type::operator[](i) = v[i];
    return *this;
  }
  
  // operator*
  array3 operator*(array3 const& x, double b) {
    array3 y(x);
    y *= b;
    return y;
  }

  array3 operator*(double b, array3 const& x) {
    return x * b;
  }

  array3 operator/(array3 const& x, double b) {
    array3 y(x);
    y /= b;
    return y;
  }

  // operator+
  array3 operator+(array3 const& x, array3 const& y) {
    array3 z(x);
    z += y;
    return z;
  }
  
  // operator-
  array3 operator-(array3 const& x, array3 const& y) {
    array3 z(x);
    z -= y;
    return z;
  }
  
  std::ostream& operator<<(std::ostream& os, array3 const& v){
    os << "( " << v[0] << " " << v[1] << " " << v[2] << " )" << std::endl;
    return os;
  }

  /* Others */
  array3 add(array3 const& x, array3 const& y, double a = 1.0) {
    array3 z( x[0] + a * y[0], x[1] + a * y[1], x[2] + a * y[2] );
    return z;
  }

  array3 cross(array3 const& x, array3 const& y) {
    array3 z( x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] );
    return z;
  }
  
} /* end namespace rpa */

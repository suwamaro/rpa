/***************************************************************************
* 3d vector (array)
*
* Copyright (C) 2018 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*
*****************************************************************************/
#ifndef _ARRAY3_H
#define _ARRAY3_H

#include <iostream>
#include <array>

namespace rpa {
  class array3 : public std::array<double, 3> {
  private:
    typedef std::array<double, 3> super_type;
  
  public:
    array3();
    array3(array3 const& v);
    explicit array3(std::initializer_list<double> init);
    explicit array3(std::size_t i);
    explicit array3(double x, double y, double z);
    void clear();
    double norm2() const;
    double norm() const;
    double operator*(array3 const& v) const;
    array3 operator*=(double b);
    array3 operator/=(double b);
    array3& operator+=(array3 const& v);
    array3& operator-=(array3 const& v);
    array3& operator=(array3 const& v);
    array3 operator-() const;    
    friend std::ostream& operator<<(std::ostream& os, array3 const& v);
  };

  /* Associated functions */
  array3 operator*(array3 const& x, double b);
  array3 operator*(double b, array3 const& x);
  array3 operator/(array3 const& x, double b);
  array3 operator+(array3 const& x, array3 const& y);
  array3 operator+(array3 const& x, array3 const& y);
  array3 operator-(array3 const& x, array3 const& y);
  std::ostream& operator<<(std::ostream& os, array3 const& v);
  array3 add(array3 const& x, array3 const& y, double a);
  array3 cross(array3 const& x, array3 const& y);
  
} /* end namespace rpa */

#endif // _ARRAY3_H

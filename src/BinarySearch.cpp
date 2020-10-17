#include <iostream>
#include "BinarySearch.h"

void BinarySearch::find_solution(double& x, double target, std::function<double(double x)> const& f, bool additive, double x_delta, double x_MIN, double x_MAX, bool debug){
  if ( debug && !std::isnan(x_MIN) && !std::isnan(x_MAX) ) { std::cout << "Finding the solution between " << x_MIN << " and " << x_MAX << std::endl; }
  
  /* Checking whether x is within the range */
  if ( ( !std::isnan( x_MIN ) && x < x_MIN ) || ( !std::isnan( x_MAX ) && x > x_MAX ) ) {
    std::cerr << "x is out of the range.\n";
    std::exit(EXIT_FAILURE);
  }

  double x2 = 0;
  if ( additive ) { x2 = x + x_delta; }
  else { x2 = 1.1 * x; }
  if ( ( !std::isnan( x_MIN ) && x2 < x_MIN ) || ( !std::isnan( x_MAX ) && x2 > x_MAX ) ) {
    std::cerr << "x2 is out of the range.\n";
    std::exit(EXIT_FAILURE);
  }

  /* Checking whether f(x) is an increasing or decreasing function of x */
  double fx = f( x );
  double fx2 = f( x2 );
  bool increasing = fx2 > fx;

  /* Determining increasing or decreasing x */
  double diff = fx - target;
  double factor = 0;
  if ( increasing ) { factor = ( diff > 0 ? 0.5 : 2. ); }
  else { factor = ( diff > 0 ? 2. : 0.5 ); }

  /* Initialization */
  // fx2 = fx;
  unsigned int n_iter = 0;
  
  /* Finding f(x) < target < f(x2) or f(x2) < target < f(x) */
  do {
    ++n_iter;

    if ( debug ) {
      std::cout << n_iter << "   " << x << "   " << fx << "   " << x2 << "  " << fx2 << "   " << target << std::endl;
    }
    
    if ( !additive && n_iter == max_iter ) {
      std::cerr << "Cannot find an optimal x.\n";
      return;
    }

    /* Storing the previous value */
    x2 = x;
    fx2 = fx;
    
    /* Updating x */
    if ( additive ) {
      if ( factor > 1. ) { x += x_delta; }
      else { x -= x_delta; }
    } else {
      x *= factor;
    }
    if ( !std::isnan( x_MIN ) && x < x_MIN ) { x = x_MIN; }
    if ( !std::isnan( x_MAX ) && x > x_MAX ) { x = x_MAX; }
    
    /* Updating the value */
    fx = f(x);

  } while( ( fx - target ) * ( fx2 - target ) > eps_fx * eps_fx );

  /* Binary search */
  if ( debug ) { std::cout << "Binary search between " << x << " and " << x2 << std::endl; }
  double dx = 0;
  if ( additive ) {
    if ( factor > 1. ) { dx = - 0.5 * x_delta; }
    else {               dx =   0.5 * x_delta; }
  } else {
    dx = 0.5 * ( 1. / factor - 1. ) * x;
  }
  diff = fx - target;

  while( std::abs( diff ) > eps_fx ) {
    ++n_iter;

    if ( debug ) {
      std::cout << n_iter << "   " << x << "   " << fx << "   " << target << "   " << diff << std::endl;
    }
    
    if ( std::abs( dx ) < 1e-24 ) {
      std::cerr << "Cannot find an optimal x.\n";
      return;
    }

    /* Updating x */
    x += dx;

    /* Updating the value */
    fx = f(x);
    diff = fx - target;
    
    /* Updating dx */
    if ( increasing ) {
      dx = ( diff > 0. ? - 0.5 * std::abs( dx ) : 0.5 * std::abs( dx ) );
    } else {
      dx = ( diff < 0. ? - 0.5 * std::abs( dx ) : 0.5 * std::abs( dx ) );
    }
    
  } /* end for while */

  if ( debug ) {
    std::cout << n_iter << "   " << x << "   " << fx << "   " << target << "   " << diff << std::endl;
  }
}

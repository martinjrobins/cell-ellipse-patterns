/*
 * GayBernePotential.cpp
 *
 *  Created on: 5 Sep 2014
 *      Author: mrobins
 */

#include "GayBernePotential.h"
#include <limits>

double GayBernePotential::evaluate(const vdouble2 &x1, const vdouble2 &u1,
                                   const vdouble2 &x2,
                                   const vdouble2 &u2) const {
  vdouble2 dx = x1 - x2;
  double r = dx.norm();
  vdouble2 rhat = dx / r;
  if (r <= std::numeric_limits<double>::min())
    rhat = vdouble2(1, 0);

  //	const double trunc = 0.000001;
  //	if (r/sigma_s < trunc) {
  //		dx = rhat*sigma_s*trunc;
  //		r = trunc*sigma_s;
  //	}

  const double udotu = u1.dot(u2);
  const double u1dotr = u1.dot(rhat);
  const double u2dotr = u2.dot(rhat);

  double q = qfunc(udotu, u1dotr, u2dotr, r);
  if (q < 0) {
    q = std::numeric_limits<double>::max();
  }
  double val =
      4 * epsilon(udotu, u1dotr, u2dotr) * (std::pow(q, 12) - std::pow(q, 6));
  if (isnan(val)) {
    return 1.0e-20 * std::numeric_limits<double>::max();
  } else {
    return val;
  }
}

double GayBernePotential::cut_off() const { return 2.5 * sigma_s * k; }

double GayBernePotential::epsilon(const double udotu, const double u1dotr,
                                  const double u2dotr) const {
  return epsilon_0 * std::pow(epsilon_dash(udotu, u1dotr, u2dotr), mu) /
         std::pow(epsilon(udotu), nu);
}

double GayBernePotential::epsilon_dash(const double udotu, const double u1dotr,
                                       const double u2dotr) const {
  const double first = std::pow(u1dotr + u2dotr, 2) / (1 + xi_dash * udotu);
  const double second = std::pow(u1dotr - u2dotr, 2) / (1 - xi_dash * udotu);
  return 1 - (xi_dash / 2.0) * (first + second);
}

double GayBernePotential::qfunc(const double udotu, const double u1dotr,
                                const double u2dotr, const double r) const {
  return sigma_s / (r - sigma(udotu, u1dotr, u2dotr) + sigma_s);
}

double GayBernePotential::epsilon(const double udotu) const {
  return sqrt(1.0 - std::pow(xi, 2) * std::pow(udotu, 2));
}

double GayBernePotential::sigma(const double udotu, const double u1dotr,
                                const double u2dotr) const {
  const double first = std::pow(u1dotr + u2dotr, 2) / (1 + xi * udotu);
  const double second = std::pow(u1dotr - u2dotr, 2) / (1 - xi * udotu);
  // if (isnan(sigma_s*sqrt(1 - (xi/2.0) * (first + second)))) std::cout <<
  // "(xi/2.0) * (first + second) = "<<(xi/2.0) * (first + second)<<" 1+xi*udotu
  // "<<1 + xi*udotu<< "1-xi*udotu "<<1 - xi*udotu <<std::endl;
  return sigma_s / sqrt(1 - (xi / 2.0) * (first + second));
}

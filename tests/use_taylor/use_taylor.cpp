/*
Copyright (c) 2009-2017 Ulf Ekstrom <uekstrom@gmail.com>

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
*/

#include <iomanip>
#include <iostream>

#include <taylor/taylor.hpp>

// This is a very small example of taking derivatives
// using libtaylor.

using namespace std;

// Define a function we want to differentiate: f(x,y)
template <typename T> T f(const T & x, const T & y) {
  return sin(log(7 * x) + exp(y)) + 9;
}

template <typename T> T slaterx(const T & x) { return -0.93 * pow(x, 4.0 / 3.0); }

int main() {
  // Compute a directional derivative of f(x,y) in the direction
  // (1,2), at point (x,y) = (3,4). This is equivalent to computing the
  // taylor expansion of f(3+1*eps, 4+2*eps) in the variable eps.
  const int Ndeg = 5; // Order of expansion.
  const int Nvar = 1; // Only one differentiation variable in this example
  taylor<double, Nvar, Ndeg> eps(0, 0); // Set up seed variable.
  taylor<double, Nvar, Ndeg> fexpansion = f(3 + eps, 4 + 2 * eps);
  // Now fexpansion contains Taylor coefficients. If we want
  // derivatives we have to multiply by the appropriate factorials
  fexpansion.deriv_facs();
  cout << "Directional derivative: " << fexpansion << endl;

  // Gradient of a 2-variable function
  taylor<double, 2, 1> dx(0.0), dy(0.0);
  dx[0] = M_E;
  dx[1] = 1.0;
  dy[0] = M_PI;
  dy[2] = 1.0;
  std::cout << dx << std::endl;
  std::cout << dy << std::endl;
  taylor<double, 2, 1> gradient = f(dx, dy);
  std::cout << std::setprecision(12) << gradient << std::endl;
  gradient.deriv_facs();
  std::cout << gradient << std::endl;

  // Hessian of a 2-variable function
  taylor<double, 2, 2> x, y;
  x[0] = M_E;
  x[1] = 1.0;
  y[0] = M_PI;
  y[2] = 1.0;
  std::cout << x << std::endl;
  std::cout << y << std::endl;
  taylor<double, 2, 2> second_order = f(x, y);
  std::cout << "Taylor expansion coefficients (grevlex)\n"
            << second_order << std::endl;
  second_order.deriv_facs();
  std::cout << "Derivatives (grevlex)\n" << second_order << std::endl;

  // Univariate case...
  taylor<double, 1, 2> dn(2.0, 0.0);
  std::cout << dn << std::endl;
  taylor<double, 1, 2> ex = slaterx(dn);
  std::cout << ex << std::endl;
  ex.deriv_facs();
  std::cout << ex << std::endl;

  return 0;
}

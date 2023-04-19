#include <boost/math/special_functions/binomial.hpp>

#include "integral.h"
#include "rys_quadrature.h"

namespace hfincpp::integral::faster_kernel::rys_quadrature {

int binomial_coefficient(const int n, int r) {
  if((r<0) || (n<r)) return 0;
  if((2*r) > n) r = n-r;
  int result = 1;
  for(int i=1; i<=r; i++) {
    result = (result * (n-i+1)) / i;
  }
  return result;
}

double binomial_expand(const int i,
                       const int a,
                       const int b,
                       const double PA,
                       const double PB) {
  double result = 0;
  for(int n=std::max(i-b, 0); n<=std::min(i, a); n++) {
    result += binomial_coefficient(a, n) * binomial_coefficient(b, i-n)
              * pow(PA, a-n) * pow(PB, b+n-i);
  }
  return result;
}


}


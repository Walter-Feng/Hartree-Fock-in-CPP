#include <boost/math/special_functions/binomial.hpp>

#include "integral.h"
#include "rys_quadrature.h"

namespace hfincpp::integral::faster_kernel::rys_quadrature {

double binomial_expand(const int i,
                       const int a,
                       const int b,
                       const double AB) {
  double result = 0;
  for(int n=std::max(i-b, 0); n<=std::min(i, a); n++) {
    result += binomial_coefficient(a, n) * binomial_coefficient(b, i-n)
              * (a == n? 1 : 0) * pow(AB, b+n-i);
  }
  return result;
}


}


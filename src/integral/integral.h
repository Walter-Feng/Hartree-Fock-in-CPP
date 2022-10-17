#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <armadillo>

#include "basis/basis.h"

namespace hfincpp::integral {

struct GaussianFunction {
  arma::vec3 center;
  arma::Col<int>::fixed<3> angular;
  double exponent;
  double coef;
};

struct ERI {
  GaussianFunction A;
  GaussianFunction B;
  GaussianFunction C;
  GaussianFunction D;
};

using GaussianFunctionPair = std::pair<GaussianFunction, GaussianFunction>;

}


#endif //INTEGRAL_H

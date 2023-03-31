#ifndef FASTER_KERNEL_INTEGRAL_H
#define FASTER_KERNEL_INTEGRAL_H

#include <armadillo>

#include "basis/basis.h"

namespace hfincpp::integral::faster_kernel {

struct UpperTriangularFunctionPairs {
  arma::vec exponents;
  arma::vec prefactors;
  arma::mat centers; // shape of n x 3
  arma::vec screening_conditions;
  arma::Mat<int> angular_momenta; // shape of n x (3 x 2)
  arma::Mat<int> indexing; // shape of n x 2, first column for i, second for j
                           // in (ij) pair
};

UpperTriangularFunctionPairs generate_function_pairs(basis::Basis basis);
}


#endif //FASTER_KERNEL_INTEGRAL_H

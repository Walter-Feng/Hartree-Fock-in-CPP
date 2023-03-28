#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <armadillo>

#include "basis/basis.h"

namespace hfincpp::integral::faster_kernel {

struct UpperTriangularShellPairs {
  arma::vec exponents;
  arma::vec prefactors;
  arma::mat centers; // shape of n x 3
  arma::vec screening_conditions;
  arma::Mat<int> angular_momenta; // shape of n x 2
  arma::Mat<int> indexing; // shape of n x 2, first column for i, second for j
                           // in (ij) pair
};

UpperTriangularShellPairs generate_shell_pairs(const basis::Basis & basis);
}


#endif //INTEGRAL_H

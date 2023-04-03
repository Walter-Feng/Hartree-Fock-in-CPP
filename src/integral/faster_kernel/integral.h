#ifndef FASTER_KERNEL_INTEGRAL_H
#define FASTER_KERNEL_INTEGRAL_H

#include <armadillo>

#include "basis/basis.h"

namespace hfincpp::integral::faster_kernel {

struct UpperTriangularShellPairs {
  arma::vec exponents; // shape of n_shell_pairs
  arma::vec prefactors;// shape of n_function_pairs
  arma::mat centers; // shape of n_shell_pairs x 3
  arma::vec screening_conditions;// shape of n_shell_pairs
  arma::Mat<int> angular_momenta; // shape of n_function_pairs x (3 x 2)
  arma::Mat<int> indexing; // shape of n_function_pairs x 2, first column for i, second for j
                           // in (ij) pair
};

UpperTriangularShellPairs generate_function_pairs(basis::Basis basis);
}


#endif //FASTER_KERNEL_INTEGRAL_H

#ifndef UTIL_ARMA_H
#define UTIL_ARMA_H

#include <armadillo>

namespace hfincpp {
namespace util {
namespace arma {

template<typename T>
::arma::Mat<T> outer_sum(const ::arma::Col<T> & col,
                         const ::arma::Col<T> & row) {
  const ::arma::uword n_cols = row.n_elem;
  const ::arma::uword n_rows = col.n_elem;
  return ::arma::repmat(col, 1, n_cols) + ::arma::repmat(row.t(), n_rows, 1);
}

}
}
}
#endif //UTIL_ARMA_H

#ifndef UTIL_GSL_CONVERTER_H
#define UTIL_GSL_CONVERTER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace gsl {

inline
gsl_vector* convert_vec(const arma::vec & vector) {
  gsl_vector * gsl_vec = gsl_vector_alloc(vector.n_elem);
  memcpy(gsl_vec->data, vector.memptr(), sizeof(double) * vector.n_elem);
  return gsl_vec;
}

inline
arma::vec convert_vec(const gsl_vector * gsl_vec) {
  arma::vec arma_vec(gsl_vec->size);
  memcpy(arma_vec.memptr(), gsl_vec->data, sizeof(double) * gsl_vec->size);
  return arma_vec;
}

inline
gsl_matrix* convert_mat(const arma::mat & mat) {
  gsl_matrix * gsl_mat = gsl_matrix_alloc(mat.n_rows, mat.n_cols);
  memcpy(gsl_mat->data, mat.memptr(),
         sizeof(double) * mat.n_rows * mat.n_cols);
  return gsl_mat;
}

inline
arma::mat convert_mat(const gsl_matrix * gsl_mat) {
  arma::mat arma_mat(gsl_mat->size1, gsl_mat->size2);
  memcpy(arma_mat.memptr(), gsl_mat->data,
         sizeof(double) * gsl_mat->size1 * gsl_mat->size2);
  return arma_mat.t();
}


}

#endif //UTIL_GSL_CONVERTER_H

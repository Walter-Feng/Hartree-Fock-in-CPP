#ifndef UTIL_GSL_CONVERTER_H
#define UTIL_GSL_CONVERTER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace gsl {

inline
gsl_vector* convert_vec(const arma::vec & vector) {

  gsl_vector * gsl_vec = gsl_vector_alloc(vector.n_elem);

#pragma omp parallel for
  for(arma::uword i=0; i<vector.n_elem; i++) {
    gsl_vector_set(gsl_vec, i, vector(i));
  }

  return gsl_vec;
}

inline
arma::vec convert_vec(const gsl_vector * gsl_vec) {

  arma::vec arma_vec(gsl_vec->size);

#pragma omp parallel for
  for(arma::uword i=0; i<arma_vec.n_elem; i++) {
    arma_vec(i) = gsl_vector_get(gsl_vec, i);
  }

  return arma_vec;
}

inline
gsl_matrix* convert_mat(const arma::mat & mat) {

  gsl_matrix * gsl_mat = gsl_matrix_alloc(mat.n_rows, mat.n_cols);

#pragma omp parallel for
  for(arma::uword i=0; i<mat.n_rows; i++) {
    for(arma::uword j=0; j<mat.n_cols; j++) {
     gsl_matrix_set(gsl_mat, i, j, mat(i,j));
    }
  }

  return gsl_mat;
}

inline
arma::mat convert_mat(const gsl_matrix * gsl_mat) {

  arma::mat arma_mat(gsl_mat->size1, gsl_mat->size2);

#pragma omp parallel for
  for(arma::uword i=0; i<arma_mat.n_rows; i++) {
    for(arma::uword j=0; j<arma_mat.n_cols; j++) {
      arma_mat(i,j) = gsl_matrix_get(gsl_mat, i, j);
    }
  }

  return arma_mat;
}


}

#endif //UTIL_GSL_CONVERTER_H

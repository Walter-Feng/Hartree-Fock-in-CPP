#include "utils.h"

namespace hfincpp::gradient::utils {

arma::cube
displacement_vectors(const arma::mat & A_xyz, const arma::mat & B_xyz) {
  arma::cube result = arma::zeros(A_xyz.n_rows, A_xyz.n_cols, B_xyz.n_cols);

  for (arma::uword i = 0; i < B_xyz.n_cols; ++i) {
    result.slice(i) = B_xyz.each_col() - B_xyz.col(i);
  }

  return result;
}

arma::cube displacement_vectors(const arma::mat & xyz) {
  arma::cube result = displacement_vectors(xyz, xyz);

  for (arma::uword i = 0; i < xyz.n_cols; ++i) {
    result.slice(i).col(i).zeros();
  }

  return result;
}

arma::cube inverse_r(const arma::mat & xyz) {
  const arma::vec norm_squared = arma::sum(arma::square(xyz)).t();
  arma::mat distance_squared = -2.0 * xyz.t() * xyz;
  distance_squared.each_col() += norm_squared;
  distance_squared.each_row() += norm_squared.t();
  arma::mat inverse_r_cube = 1.0 / arma::pow(distance_squared, 1.5);
  inverse_r_cube.diag().ones();

  const arma::cube displacement = displacement_vectors(xyz);

  arma::cube result(xyz.n_cols, xyz.n_cols, xyz.n_rows);
  for (arma::uword i = 0; i < xyz.n_rows; i++) {
    const arma::mat displacement_in_this_dimension = displacement.row(i);
    result.slice(i) = displacement_in_this_dimension % inverse_r_cube;
  }

  return result;

}
}
#include <boost/math/special_functions/binomial.hpp>

#include "integral.h"

#include "util/armadillo.h"

namespace hfincpp::integral::faster_kernel {

template<typename T>
arma::Mat<T> concatenate_vertically(std::vector<arma::Mat<T>> & matrices) {
  arma::Mat<T> head = matrices[0];
  for(size_t i=0; i<matrices.size(); i++) {
    head = arma::join_vert(head, matrices[i]);
  }

  return head;
}

UpperTriangularShellPairs generate_shell_pairs(const basis::Basis & basis) {

  arma::uword n_shell_pairs = basis.n_functions() * (basis.n_functions() + 1) / 2;
  std::vector<UpperTriangularShellPairs> list_of_shell_pairs[n_shell_pairs];

#pragma omp parallel for collapse(2)
  for (int i = 0; i < basis.n_functions(); i++) {
    for (int j = i; j < basis.n_functions(); j++) {
      const auto & function_i = basis.functions[i];
      const auto n_gto_from_i = function_i.coefficients.n_elem;
      const auto & function_j = basis.functions[j];
      const auto n_gto_from_j = function_j.coefficients.n_elem;

      const size_t upper_triangular_indexing =
          (basis.n_functions() + basis.n_functions() - i + 1) * i / 2;

      const int n_binomial_terms =
          arma::prod(function_i.angular + function_j.angular + 1);

      const arma::vec new_exponents =
          arma::repmat(
              arma::vectorise(
                  util::arma::outer_sum(function_i.exponents,
                                        function_j.exponents)),
          n_binomial_terms, 1);

      const double center_distance_squared =
          arma::sum(arma::square(function_i.center - function_j.center));

      const arma::vec multiplied_exponents =
          arma::repmat(
              arma::vectorise(function_i.exponents * function_j.exponents.t())
          , n_binomial_terms, 1);

      const arma::vec exponential_prefactor =
          arma::exp(- center_distance_squared * multiplied_exponents / new_exponents);

      arma::mat new_center =
          function_i.exponents * function_i.center.t() +
              function_j.exponents * function_j.center.t();

      new_center.each_col() /= new_exponents;



      arma::vec3 binomial_coefficients;


    }
  }
}

}
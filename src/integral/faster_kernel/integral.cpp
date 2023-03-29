#include <boost/math/special_functions/binomial.hpp>

#include "integral.h"

#include "integral/rys_quadrature.h"
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

UpperTriangularFunctionPairs generate_function_pairs(const basis::Basis & basis) {

  arma::uword n_function_pairs = basis.n_functions() * (basis.n_functions() + 1) / 2;
  std::vector<UpperTriangularFunctionPairs> list_of_function_pairs(n_function_pairs);

#pragma omp parallel for collapse(2)
  for (int i = 0; i < basis.n_functions(); i++) {
    for (int j = i; j < basis.n_functions(); j++) {
      const auto & function_i = basis.functions[i];
      const auto n_gto_from_i = function_i.coefficients.n_elem;
      const auto & function_j = basis.functions[j];
      const auto n_gto_from_j = function_j.coefficients.n_elem;

      const auto n_gto = n_gto_from_i * n_gto_from_j;

      const size_t upper_triangular_indexing =
          (basis.n_functions() + basis.n_functions() - i + 1) * i / 2 + j - i;

      const arma::vec new_exponents =
              arma::vectorise(
                  util::arma::outer_sum(function_i.exponents,
                                           function_j.exponents));

      const double center_distance_squared =
          arma::sum(arma::square(function_i.center - function_j.center));

      const arma::vec multiplied_exponents =
          arma::vectorise(function_i.exponents * function_j.exponents.t());

      const arma::vec exponential_prefactor =
          arma::exp(- center_distance_squared * multiplied_exponents / new_exponents);

      arma::mat new_center(n_gto_from_i * n_gto_from_j, 3);

      const arma::mat i_exponent_weighted_center =
          function_i.exponents * function_i.center.t();
      const arma::mat j_exponent_weighted_center =
          function_j.exponents * function_j.center.t();

      for(int xyz_index=0; xyz_index<3; xyz_index++) {
        const arma::vec i_share = i_exponent_weighted_center.col(xyz_index);
        const arma::vec j_share = j_exponent_weighted_center.col(xyz_index);

        const arma::mat outer_sum = util::arma::outer_sum(i_share, j_share);
        new_center.col(xyz_index) =
               arma::vectorise(outer_sum) / new_exponents;
      }

      const arma::vec new_coefficients =
          arma::vectorise(function_i.coefficients * function_j.coefficients.t());

      const arma::Mat<int> angular_momenta =
          arma::repmat(arma::join_vert(function_i.angular,
                                       function_j.angular).t(), n_gto, 1);

      const arma::Mat<int> indexing =
          arma::repmat(arma::Mat<int>{{i, j}}, n_gto, 1);

      arma::mat screening_condition(n_gto_from_i, n_gto_from_j);
      for(arma::uword gto_i = 0; gto_i < n_gto_from_i; gto_i++) {
        for(arma::uword gto_j = 0; gto_j < n_gto_from_j; gto_j++) {
          const GaussianFunction gto_function_i{function_i.center,
                                                function_i.angular,
                                                function_i.exponents(
                                                    gto_i),
                                                function_i.coefficients(
                                                    gto_i)};

          const GaussianFunction gto_function_j{function_j.center,
                                                function_j.angular,
                                                function_j.exponents(
                                                    gto_j),
                                                function_j.coefficients(
                                                    gto_j)};

          const ERI eri_info{gto_function_i, gto_function_j,
                             gto_function_i, gto_function_j};

          screening_condition(gto_i, gto_j) =
              std::sqrt(rys_quadrature::electron_repulsive_integral(eri_info));
        }
      }

      const arma::vec prefactors = exponential_prefactor % new_coefficients;

      list_of_function_pairs[upper_triangular_indexing] =
          {new_exponents,
           prefactors,
           new_center,
           arma::vectorise(screening_condition) % prefactors,
           angular_momenta,
           indexing};
    }
  }

  UpperTriangularFunctionPairs result;
  for(const auto & i_pair : list_of_function_pairs) {
    result.exponents = arma::join_vert(result.exponents, i_pair.exponents);
    result.prefactors = arma::join_vert(result.prefactors, i_pair.prefactors);
    result.centers = arma::join_vert(result.centers, i_pair.centers);
    result.screening_conditions = arma::join_vert(result.screening_conditions,
                                                  i_pair.screening_conditions);
    result.angular_momenta = arma::join_vert(result.angular_momenta,
                                             i_pair.angular_momenta);
    result.indexing = arma::join_vert(result.indexing, i_pair.indexing);
  }

  return result;
}

}
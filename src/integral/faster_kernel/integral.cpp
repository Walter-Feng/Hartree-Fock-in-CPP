#include "integral.h"

#include "integral/rys_quadrature.h"
#include "util/armadillo.h"

namespace hfincpp::integral::faster_kernel {

UpperTriangularShellPairs generate_shell_pairs(basis::Basis basis) {

  basis = basis.sort_by_angular_momentum();
  arma::uword n_shell_pairs = basis.n_shells() * (basis.n_shells() + 1) / 2;
  std::vector<UpperTriangularShellPairs> list_of_function_pairs(n_shell_pairs);

  std::vector<arma::uword> n_functions_per_shell;

#pragma omp parallel for collapse(2)
  for (int i = 0; i < basis.n_shells(); i++) {
    for (int j = i; j < basis.n_shells(); j++) {
      const auto & shell_i = basis.shells[i];
      const auto n_gto_from_i = shell_i.exponents.n_elem;
      const auto & shell_j = basis.shells[j];
      const auto n_gto_from_j = shell_j.exponents.n_elem;

      const auto n_gto = n_gto_from_i * n_gto_from_j;

      const size_t upper_triangular_indexing =
          (basis.n_functions() + basis.n_functions() - i + 1) * i / 2 + j - i;

      const arma::vec new_exponents =
              arma::vectorise(
                  util::arma::outer_sum(shell_i.exponents,
                                        shell_j.exponents));

      const double center_distance_squared =
          arma::sum(arma::square(shell_i.center - shell_j.center));

      const arma::vec multiplied_exponents =
          arma::vectorise(shell_i.exponents * shell_j.exponents.t());

      const arma::vec exponential_prefactor =
          arma::exp(- center_distance_squared * multiplied_exponents / new_exponents);

      arma::mat new_center(n_gto_from_i * n_gto_from_j, 3);

      const arma::mat i_exponent_weighted_center =
          shell_i.exponents * shell_i.center.t();
      const arma::mat j_exponent_weighted_center =
          shell_j.exponents * shell_j.center.t();

      for(int xyz_index=0; xyz_index<3; xyz_index++) {
        const arma::vec i_share = i_exponent_weighted_center.col(xyz_index);
        const arma::vec j_share = j_exponent_weighted_center.col(xyz_index);

        const arma::mat outer_sum = util::arma::outer_sum(i_share, j_share);
        new_center.col(xyz_index) =
               arma::vectorise(outer_sum) / new_exponents;
      }

      const arma::uvec function_indices_for_shell_i =
          arma::find(basis.shell_indices == i);
      const arma::uvec function_indices_for_shell_j =
          arma::find(basis.shell_indices == j);

      n_functions_per_shell.push_back(
          function_indices_for_shell_i.n_elem * function_indices_for_shell_j.n_elem);

      arma::vec new_coefficients;
      arma::Mat<int> angular_momenta;
      double screening_condition;

      for(arma::uword function_j=0;
          function_j<function_indices_for_shell_j.n_elem;
          function_j++) {
        for(arma::uword function_i=0;
            function_i<function_indices_for_shell_i.n_elem;
            function_i++) {
          const auto & i_function_object = basis.functions[function_i];
          const auto & j_function_object = basis.functions[function_j];
          new_coefficients =
              arma::join_vert(new_coefficients,
                              arma::vectorise(
                                  i_function_object.coefficients
                                  * j_function_object.coefficients.t()));

          angular_momenta = arma::join_vert(
              angular_momenta,
              arma::repmat(arma::join_vert(i_function_object.angular,
                                           j_function_object.angular).t(), n_gto, 1)
              );


        }
      }


      const arma::Mat<int> indexing =
          arma::repmat(arma::Mat<int>{{i, j}}, n_gto, 1);

      for(arma::uword gto_i = 0; gto_i < n_gto_from_i; gto_i++) {
        for(arma::uword gto_j = 0; gto_j < n_gto_from_j; gto_j++) {
          const GaussianFunction gto_function_i{shell_i.center,
                                                shell_i.angular,
                                                shell_i.exponents(
                                                    gto_i),
                                                shell_i.coefficients(
                                                    gto_i)};

          const GaussianFunction gto_function_j{shell_j.center,
                                                shell_j.angular,
                                                shell_j.exponents(
                                                    gto_j),
                                                shell_j.coefficients(
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

  UpperTriangularShellPairs result;
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
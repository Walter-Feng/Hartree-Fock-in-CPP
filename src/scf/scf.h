#ifndef HFINCPP_SCF_H
#define HFINCPP_SCF_H

#include <armadillo>
#include <functional>
#include <json.hpp>

#include "util/json.h"
#include "util/time.h"
#include "util/printer.h"
#include "global/error.h"

namespace hfincpp::scf {

template<typename T>
using FockMatrix = arma::Cube<T>;

template<typename T>
using OverlapMatrix = arma::Cube<T>;

template<typename T>
using DensityMatrix = arma::Cube<T>;

using Eigenvalue = arma::vec;
using OccupationVector = arma::vec;

template<typename T>
using EnergyBuilder = std::function<double(const DensityMatrix<T> &)>;

using OccupationBuilder = std::function<OccupationVector(const arma::vec &,
                                                         const double)>;

template<typename T>
using FockBuilder = std::function<FockMatrix<T>(const DensityMatrix<T> &)>;

struct SimpleSCFWrapper {
  double energy;
  double diff;
};

template<typename SCFWrapper>
Printer<SCFWrapper> generic_scf_printer = [](const SCFWrapper & state,
                                             const double computation_time,
                                             const int iter,
                                             const int print_level,
                                             const bool print_header) -> int {

  int width = 18;
  int precision = 8;

  if (print_level >= 2) {
    width = 27;
    precision = 17;
  }

  int total_length = 0;

  if (print_level >= 1) {
    total_length = 6 + width * 3;

    if (print_header) {
      fmt::print("Performing SCF ...\n");
      print_separator(total_length);
      fmt::print("{:>{}}", "|Iter|", 6);
      fmt::print("{:>{}}", "Time / s |", width);
      fmt::print("{:>{}}", "Energy / a.u. |", width);
      fmt::print("{:>{}}", "Energy Diff |", width);
      fmt::print("\n");
      print_separator(total_length);
    }

    fmt::print("{:>{}}", iter, 6);
    const arma::rowvec combined{computation_time, state.energy, state.diff};
    print(combined, width, precision);
    fmt::print("\n");
  }

  return total_length;

};


template<typename T>
struct SCFResult {
  arma::mat eigenvalues;
  arma::Cube<T> orbitals;
  arma::mat occupations;
  DensityMatrix<T> density;
  OverlapMatrix<T> overlap;
  FockMatrix<T> fock;

  double energy;

  [[nodiscard]] nlohmann::json to_json() const {

    const auto n_items = eigenvalues.n_cols;
    assert(orbitals.n_slices == n_items);
    assert(occupations.n_cols == n_items);
    assert(density.n_slices == n_items);
    assert(overlap.n_slices == n_items);
    assert(fock.n_slices == n_items);

    if (n_items == 1) {
      nlohmann::json result;
      util::put(result, "energy", energy);
      util::put(result, "eigenvalues", eigenvalues.col(0));
      util::put(result, "orbitals", orbitals.slice(0));
      util::put(result, "occupations", occupations.col(0));
      util::put(result, "density", density.slice(0));
      util::put(result, "overlap", overlap.slice(0));
      util::put(result, "fock", fock.slice(0));

      return result;
    } else {
      nlohmann::json::array_t array;
      for (arma::uword i = 0; i < n_items; i++) {
        nlohmann::json channel;
        util::put(channel, "eigenvalues", eigenvalues.col(0));
        util::put(channel, "orbitals", orbitals.slice(0));
        util::put(channel, "occupations", occupations.col(0));
        util::put(channel, "density", density.slice(0));
        util::put(channel, "overlap", overlap.slice(0));
        util::put(channel, "fock", fock.slice(0));

        array.push_back(array);
      }

      nlohmann::json result;
      result["channels"] = array;
      result["energy"] = energy;

      return result;
    }

    __builtin_unreachable();
  }

};

template<class State>
using Update = std::function<State(const State & state)>;

template<class UpdateMethod, typename T>
SCFResult<T> scf(const EnergyBuilder<T> & energy_builder,
                 const FockBuilder<T> & fock_builder,
                 const OccupationBuilder & occupation_builder,
                 const UpdateMethod & update_method,
                 const OverlapMatrix<T> & overlap,
                 const DensityMatrix<T> & initial,
                 const arma::vec & n_electrons,
                 const int max_iter,
                 const double energy_tolerance,
                 const int print_level = 1) {

  const double initial_energy = energy_builder(initial);

  SimpleSCFWrapper wrapper = {initial_energy, 0};
  const int total_length =
      generic_scf_printer<SimpleSCFWrapper>(wrapper, 0.0, 0, print_level, true);

 Timer scf_time;

  UpdateMethod updater = update_method;
  DensityMatrix<T> previous_state = initial;
  double previous_energy = initial_energy;

  int iter = 1;

  const auto n_ao = overlap.slice(0).n_rows;

  for (; iter <= max_iter; iter++) {
    std::pair<Update<DensityMatrix<T>>, UpdateMethod> renewed = updater(
        previous_state);

    const DensityMatrix<T> updated_state = renewed.first(previous_state);
    updater = renewed.second;
    const FockMatrix<T> fock_matrices = fock_builder(updated_state);

    DensityMatrix<T> new_state(n_ao, n_ao, fock_matrices.n_slices);
    arma::mat eigenvalues(n_ao, fock_matrices.n_slices);
    arma::Cube<T> orbitals(n_ao, n_ao, fock_matrices.n_slices);
    arma::mat occupation_vectors(n_ao, fock_matrices.n_slices);
    for (arma::uword i = 0; i < fock_matrices.n_slices; i++) {
      arma::cx_vec eigvals;
      arma::cx_mat eigvecs;

      arma::eig_pair(eigvals, eigvecs, fock_matrices.slice(i),
                     overlap.slice(i));

      const arma::vec real_eigenvalues = arma::real(eigvals);
      const arma::uvec sort_index = arma::sort_index(real_eigenvalues);
      const arma::vec sorted_eigenvalues = real_eigenvalues(sort_index);
      arma::cx_mat sorted_eigvecs = eigvecs.cols(sort_index);
      const arma::rowvec normalization_constant =
          arma::real(arma::sum(sorted_eigvecs % (overlap.slice(0) * sorted_eigvecs)));

      sorted_eigvecs.each_row() %=
          1.0 / arma::conv_to<arma::cx_mat>::from(arma::sqrt(normalization_constant));

      eigenvalues.col(i) = sorted_eigenvalues;

      if constexpr(std::is_same<T, double>::value) {
        orbitals.slice(i) = arma::real(sorted_eigvecs);
      } else {
        orbitals.slice(i) = sorted_eigvecs;
      }

      const OccupationVector occupation =
          occupation_builder(sorted_eigenvalues, n_electrons(i));

      occupation_vectors.col(i) = occupation;

      const arma::Mat<T> new_density =
          orbitals.slice(i) * arma::diagmat(occupation) * orbitals.slice(i).t();

      new_state.slice(i) = new_density;
    }
    double new_energy = energy_builder(new_state);
    double diff = new_energy - previous_energy;
    wrapper = {new_energy, diff};
    generic_scf_printer<SimpleSCFWrapper>(wrapper, scf_time.elapsed(), iter,
                                          print_level, false);
    if (std::abs(diff) < energy_tolerance) {
      if(print_level >= 1) {
        print_separator(total_length);
      }

      return {eigenvalues, orbitals, occupation_vectors, new_state, overlap,
              fock_matrices, new_energy};
    }
    previous_state = new_state;
    previous_energy = new_energy;
  }

  throw Error("SCF did not converge.");
}

}

#endif //HFINCPP_SCF_H

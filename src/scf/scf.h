#ifndef HFINCPP_SCF_H
#define HFINCPP_SCF_H

#include <armadillo>
#include <functional>
#include <chrono>
#include <json.hpp>

#include "util/printer.h"
#include "util/error.h"

namespace hfincpp::scf {

template<typename T>
using FockMatrix = std::vector<arma::Mat<T>>;

template<typename T>
using OverlapMatrix = std::vector<arma::Mat<T>>;

template<typename T>
using DensityMatrix = std::vector<arma::Mat<T>>;

using Eigenvalue = arma::vec;
using OccupationVector = arma::vec;

template<typename T>
using EnergyBuilder = std::function<double(const DensityMatrix<T> &)>;

using OccupationBuilder = std::function<OccupationVector(const arma::vec &,
                                                         const double)>;

template<typename T>
using FockBuilder = std::function<FockMatrix<T>(const DensityMatrix<T> &)>;

template<typename T>
struct SCFWrapper {
  double energy;
  double diff;
};

template<typename T>
struct SCFResult{
  std::vector<Eigenvalue> eigenvalues;
  std::vector<arma::Mat<T>> orbitals;
  std::vector<OccupationVector> occupations;
  DensityMatrix<T> density;
  OverlapMatrix<T> overlap;
  FockMatrix<T> fock;

  double energy;
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

  SCFWrapper<T> wrapper = {initial, initial_energy, 0};
  const int total_length = generic_scf_printer(wrapper, 0, 0, print_level, true);

  auto start_scf = std::chrono::high_resolution_clock::now();

  UpdateMethod updater = update_method;
  DensityMatrix<T> previous_state = initial;
  double previous_energy = initial_energy;

  int iter = 1;
  for(;iter <= max_iter; iter++) {
    std::pair<DensityMatrix<T>, UpdateMethod> renewed = updater(previous_state);
    const DensityMatrix<T> updated_state = renewed.first(previous_state);
    updater = renewed.second;
    const FockMatrix<T> fock_matrices = fock_builder(updated_state);

    DensityMatrix<T> new_state;
    std::vector<Eigenvalue> eigenvalues;
    std::vector<arma::Mat<T>> orbitals;
    std::vector<OccupationVector> occupation_vectors;
    for(int i=0; i<fock_matrices.size(); i++) {
      arma::Col<T> eigvals;
      arma::Mat<T> eigvecs;

      arma::eig_pair(eigvals, eigvecs, fock_matrices[i], overlap[i]);

      const arma::vec real_eigenvalues = arma::real(eigvals);

      eigenvalues.push_back(real_eigenvalues);
      orbitals.push_back(eigvecs);

      const OccupationVector occupation =
          occupation_builder(real_eigenvalues, n_electrons(i));

      occupation_vectors.push_back(occupation);

      const arma::Mat<T> new_density =
          eigvecs * arma::diagmat(occupation) * eigvecs.t();

      new_state.push_back(new_density);
    }
    const double new_energy = energy_builder(new_state);
    const double diff = new_energy - previous_energy;
    wrapper = {new_state, new_energy, diff};
    auto now = std::chrono::high_resolution_clock::now();
    generic_scf_printer(wrapper, iter, now - start_scf, print_level);
    if(std::abs(diff) < energy_tolerance) {

      for (int i = 0; i < total_length; i++) {
        fmt::print("=");
      }
      fmt::print("\n");

      fmt::print("\n");

      return {eigenvalues, orbitals, occupation_vectors, new_state, overlap, fock_matrices, new_energy};
    }
    previous_state = new_state;
    previous_energy = new_energy;
  }

  throw Error("SCF did not converge.");
}


template<typename SCFWrapper, typename T>
Printer<SCFWrapper> generic_scf_printer = [](const SCFWrapper & state,
                                        const double computation_time,
                                        const int iter,
                                        const int print_level = 1,
                                        const bool print_header = false) -> int {

  int width = 18;
  int precision = 8;

  if (print_level > 2) {
    width = 27;
    precision = 17;
  }

  int total_length = 0;

  if (print_level == 1) {
    const double energy = state.positional_expectation();
    total_length = 6 + width * 2;

    if (print_header) {
      for (int i = 0; i < total_length; i++) {
        fmt::print("=");
      }
      fmt::print("\n");

      fmt::print("{:>{}}", "|Iter|", 6);
      fmt::print("{:>{}}", "Time |", width);
      fmt::print("{:>{}}", "Energy |", width);
      fmt::print("{:>{}}", "Energy Diff |", width);
      fmt::print("\n");
      for (int i = 0; i < total_length; i++) {
        fmt::print("=");
      }
      fmt::print("\n");
    }

    fmt::print("{:>{}}", iter, 6);
    const arma::rowvec combined{computation_time, state.energy, state.diff};
    print(combined, width, precision);
  }

  return total_length;

};

}

#endif //HFINCPP_SCF_H

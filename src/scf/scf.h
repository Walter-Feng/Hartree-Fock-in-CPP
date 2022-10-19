#ifndef HFINCPP_SCF_H
#define HFINCPP_SCF_H

#include <armadillo>
#include <functional>
#include <chrono>
#include <json.hpp>

#include "util/json.h"
#include "util/printer.h"
#include "util/error.h"

namespace hfincpp::scf {

template<typename T>
using FockMatrix = std::vector<arma::Mat<T>>;

template<typename T>
using OverlapMatrix = std::vector<arma::Mat<T>>;

template<typename T>
using DensityMatrix = std::vector<arma::Mat<T>>;

template<typename T>
std::vector<arma::Mat<T>> operator+(const std::vector<arma::Mat<T>> & A,
                                     const std::vector<arma::Mat<T>> & B) {
  assert(A.size() == B.size());

  std::vector<arma::Mat<T>> result(A.size());

#pragma omp parallel for
  for(unsigned long i=0; i<A.size(); i++) {
    result[i] = A[i] + B[i];
  }

  return result;
}

template<typename T>
std::vector<arma::Mat<T>> operator*(const std::vector<arma::Mat<T>> & A,
                                    const double factor) {
  std::vector<arma::Mat<T>> result(A.size());

#pragma omp parallel for
  for(unsigned long i=0; i<A.size(); i++) {
    result[i] = A[i] * factor;
  }

  return result;
}

template<typename T>
std::vector<arma::Mat<T>> operator*(const double factor,
                                    const std::vector<arma::Mat<T>> & A) {
  std::vector<arma::Mat<T>> result(A.size());

#pragma omp parallel for
  for(unsigned long i=0; i<A.size(); i++) {
    result[i] = A[i] * factor;
  }

  return result;
}

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

  if (print_level > 2) {
    width = 27;
    precision = 17;
  }

  int total_length = 0;

  if (print_level == 1) {
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


template<typename T>
struct SCFResult{
  std::vector<Eigenvalue> eigenvalues;
  std::vector<arma::Mat<T>> orbitals;
  std::vector<OccupationVector> occupations;
  DensityMatrix<T> density;
  OverlapMatrix<T> overlap;
  FockMatrix<T> fock;

  double energy;

  [[nodiscard]] nlohmann::json to_json() const {

    const auto n_items = eigenvalues.size();
    assert(orbitals.size() == n_items);
    assert(occupations.size() == n_items);
    assert(density.size() == n_items);
    assert(overlap.size() == n_items);
    assert(fock.size() == n_items);

    if(n_items == 1) {
      nlohmann::json result;
      util::put(result, "eigenvalues", eigenvalues[0]);
      util::put(result, "orbitals", orbitals[0]);
      util::put(result, "occupations", occupations[0]);
      util::put(result, "density", density[0]);
      util::put(result, "overlap", overlap[0]);
      util::put(result, "fock", fock[0]);

      return result;
    } else {
      nlohmann::json::array_t array;
      for(int i=0; i<n_items; i++) {
        nlohmann::json channel;
        util::put(channel, "eigenvalues", eigenvalues[i]);
        util::put(channel, "orbitals", orbitals[i]);
        util::put(channel, "occupations", occupations[i]);
        util::put(channel, "density", density[i]);
        util::put(channel, "overlap", overlap[i]);
        util::put(channel, "fock", fock[i]);

        array.push_back(array);

        return array;
      }
    }
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

  auto start_scf = std::chrono::high_resolution_clock::now();

  UpdateMethod updater = update_method;
  DensityMatrix<T> previous_state = initial;
  double previous_energy = initial_energy;

  int iter = 1;
  for(;iter <= max_iter; iter++) {
    std::pair<Update<DensityMatrix<T>>, UpdateMethod> renewed = updater(previous_state);
    const DensityMatrix<T> updated_state = renewed.first(previous_state);
    updater = renewed.second;
    const FockMatrix<T> fock_matrices = fock_builder(updated_state);

    DensityMatrix<T> new_state;
    std::vector<Eigenvalue> eigenvalues;
    std::vector<arma::Mat<T>> orbitals;
    std::vector<OccupationVector> occupation_vectors;
    for(int i=0; i<fock_matrices.size(); i++) {
      arma::cx_vec eigvals;
      arma::cx_mat eigvecs;

      arma::eig_pair(eigvals, eigvecs, fock_matrices[i], overlap[i]);

      const arma::vec real_eigenvalues = arma::real(eigvals);

      eigenvalues.push_back(real_eigenvalues);
      if constexpr(std::is_same<T, double>::value) {
        orbitals.push_back(arma::real(eigvecs));
      } else {
        orbitals.push_back(eigvecs);
      }


      const OccupationVector occupation =
          occupation_builder(real_eigenvalues, n_electrons(i));

      occupation_vectors.push_back(occupation);

      const arma::Mat<T> new_density =
          eigvecs * arma::diagmat(occupation) * eigvecs.t();

      new_state.push_back(new_density);
    }
    double new_energy = energy_builder(new_state);
    double diff = new_energy - previous_energy;
    wrapper = {new_energy, diff};
    auto now = std::chrono::high_resolution_clock::now();
    const double time_consumed = (now - start_scf).count();
    generic_scf_printer<SimpleSCFWrapper>(wrapper, time_consumed, iter, print_level, false);
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

}

#endif //HFINCPP_SCF_H

#ifndef HFINCPP_BASIS_H
#define HFINCPP_BASIS_H

#include <armadillo>

#include <string>

#include "geometry/geometry.h"

namespace hfincpp {
namespace basis {

struct GTOFunction {
  arma::vec3 center;
  arma::Col<int>::fixed<3> angular;
  arma::vec exponents;
  arma::vec coefficients;
};

struct GTOShell {
  arma::vec3 center;
  int angular;
  arma::vec exponents;
};


struct Basis {
  std::vector<GTOFunction> functions;
  std::vector<GTOShell> shells;
  std::vector<std::string> atom_symbols;
  arma::uvec atomic_numbers;
  arma::uvec atom_indices;
  arma::uvec shell_indices;
  std::vector<std::string> function_labels;
  std::vector<std::string> shell_labels;
  std::string basis_name;
  arma::uvec to_conventional_function_indexing;

  Basis();

  Basis(const geometry::Atoms & atoms,
        const std::string & basis_name);

  int n_atoms() const;
  int n_shells() const;
  int n_functions() const;

  [[nodiscard]] arma::uvec on_atom(const arma::uword atom_index) const;
  std::vector<arma::uvec> on_atoms() const;
  Basis sort_by_angular_momentum() const;
};
}
}


#endif //HFINCPP_BASIS_H

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

struct Basis {
  std::vector<GTOFunction> functions;
  std::vector<std::string> atom_symbols;
  arma::uvec atomic_numbers;
  std::vector<std::string> function_labels;

  Basis();

  Basis(const geometry::Atoms & atoms,
        const std::string & basis_name);

  Basis(const Basis & basis);

  int n_atoms() const;
  int n_functions() const;

};
}
}


#endif //HFINCPP_BASIS_H

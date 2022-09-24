#ifndef HFINCPP_BASIS_H
#define HFINCPP_BASIS_H

#include <armadillo>

#include <string>

namespace hfincpp {
namespace basis {

struct GTOFunction {
  arma::vec3 center;
  arma::Col<int>::fixed<3> angular;
  arma::vec exponents;
  arma::vec coefs;
};

struct GTOShell {
  arma::vec3 center;
  std::vector<std::pair<arma::Col<int>::fixed<3>, double>> angular_momentum;
  arma::vec exponents;
  arma::vec coefs;
  std::string atom_symbol;
  int atom_number;

  GTOShell gradient(const int xyz_index);
};

struct Basis {
  std::vector<GTOShell> shells;
  std::vector<std::string> atom_symbols;
  arma::uvec atom_numbers;

  int n_atoms() const;
  int n_shells() const;
  std::vector<std::string> labels() const;


};
}
}


#endif //HFINCPP_BASIS_H

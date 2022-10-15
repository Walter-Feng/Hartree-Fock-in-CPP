#include "basis.h"

#include "geometry/periodic_table.h"
#include "util/json.h"
#include "util/error.h"

namespace hfincpp {
namespace basis {

GTOShell GTOShell::gradient(const int xyz_index) {

  std::vector<std::pair<arma::Col<int>::fixed<3>, double>> new_angular_momentum;
  for(int i=0; i<angular_momentum.size(); i++) {
    if(angular_momentum[i].first(xyz_index) > 0) {
      arma::Col<int>::fixed<3> derived_angular_momentum = angular_momentum[i].first;
      derived_angular_momentum(xyz_index) -= 1;

      new_angular_momentum.push_back(
          {derived_angular_momentum, xyz_index * derived_angular_momentum(xyz_index)});
    }
  }

  return {center, new_angular_momentum, exponents, coefs, atom_symbol, atom_number};
}

Basis::Basis() {shells = {}; atom_symbols = {}; atom_numbers = arma::uvec{};}

Basis::Basis(const geometry::Atoms & atoms,
             const std::string & basis_name) {
  std::ifstream f("../data/" + basis_name + ".0.json");
  if(f.is_open()) {
    nlohmann::json data = nlohmann::json::parse(f);

    Basis basis;

    const int n_atoms = atoms.n_atoms();
    for(int i=0; i<n_atoms; i++) {
      const int atomic_numbers = atoms.atomic_numbers[i];
      
    }

  } else {
    throw Error("The basis is not found");
  }
}

Basis::Basis(const Basis & basis) {
shells = basis.shells;
atom_symbols = basis.atom_symbols;
atom_numbers = basis.atom_numbers;
}

}
}
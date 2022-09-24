#include "basis.h"

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
}
}
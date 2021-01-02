#include "basis.h"

namespace hfincpp {
namespace basis {

GTOShell GTOShell::gradient(const int xyz_index) {

  std::vector<std::pair<arma::ivec3, double>> new_angular_momentum;
  for(int i=0; i<angular_momentum.size(); i++) {
    if(angular_momentum[i].first(xyz_index) > 0) {
      arma::ivec3 derived_angular_momentum = angular_momentum[i].first;
      derived_angular_momentum(xyz_index) -= 1;

      new_angular_momentum.push_back(
          {derived_angular_momentum, xyz_index * derived_angular_momentum(xyz_index)});
    }
  }

  return {center, new_angular_momentum, GTO_alpha, GTO_coef, atom_symbol, atom_number};
}
}
}
#ifndef HFINCPP_GEOMETRY_H
#define HFINCPP_GEOMETRY_H

#include <string>
#include <armadillo>
#include <vector>

namespace hfincpp {

namespace geometry {
struct Atoms {

  std::vector<std::string> symbols;
  arma::uvec atomic_numbers;
  arma::mat xyz;
  int charge = 0;

  int n_atoms() const;
  int n_elec() const;

};

}
}

#endif //HFINCPP_GEOMETRY_H

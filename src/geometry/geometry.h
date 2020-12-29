#ifndef HFINCPP_GEOMETRY_H
#define HFINCPP_GEOMETRY_H

#include <string>
#include <armadillo>
#include <vector>

namespace hfincpp {

namespace geometry {
  struct atoms {

    std::vector<std::string> symbols;
    arma::uvec atomic_numbers;
    arma::mat xyz;

  };

  struct geometry {
    atoms atoms;
    int charge;

    //(TODO) Point Group
  };

}
}

#endif //HFINCPP_GEOMETRY_H

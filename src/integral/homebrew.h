#ifndef HFINCPP_HOMEBREW_INTEGRAL_H
#define HFINCPP_HOMEBREW_INTEGRAL_H

#include <armadillo>
#include <basis/basis.h>

namespace hfincpp {
namespace integral {

class HomebrewEngine {

private:
  arma::vec two_electron_tensor;
  arma::vec overlap;
  arma::vec kinetic;
  arma::vec nuclear_attraction;

public:

  HomebrewEngine(const basis::Basis & basis);

  double overlap_integral(const int shell_index_i,
                          const int shell_index_j) const;

  double kinetic_integral(const int shell_index_i,
                          const int shell_index_j) const;

  double nuclear_attraction_integral(const int shell_index_i,
                                     const int shell_index_j) const;

  double coulomb_integral(const int shell_index_i,
                          const int shell_index_j,
                          const int shell_index_k,
                          const int shell_index_l) const;
};
}
}





#endif //HFINCPP_HOMEBREW_INTEGRAL_H

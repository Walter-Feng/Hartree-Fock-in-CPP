#include "geometry.h"

#include <cassert>

namespace hfincpp::geometry {

int Atoms::n_atoms() const {
  assert(atomic_numbers.n_elem == xyz.n_cols);
  assert(atomic_numbers.n_elem == symbols.size());

  return atomic_numbers.n_elem;
}

int Atoms::n_elec() const {
  return arma::sum(atomic_numbers) - charge;
}

}
#include "occupation.h"

namespace hfincpp::scf::occupation {

OccupationBuilder simple_occupation(const double n_elec_per_orb) {
  return [n_elec_per_orb](const arma::vec & eigenvalues,
                          const double n_elec) {
    assert(eigenvalues.is_sorted());

    arma::vec occupation_vector(arma::size(eigenvalues), arma::fill::zeros);
    for (int i = 0; i < std::floor(n_elec / n_elec_per_orb); i++) {
      occupation_vector(i) = n_elec_per_orb;
    }
    occupation_vector(std::ceil(n_elec / n_elec_per_orb)) =
        n_elec - arma::sum(occupation_vector);

    return occupation_vector;
  };
}

}
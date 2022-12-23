#ifndef HF_RHF_H
#define HF_RHF_H

#include <json.hpp>
#include "geometry/geometry.h"
#include "basis/basis.h"
#include "scf/scf.h"

namespace hfincpp::hf {

arma::mat core_hamiltonian(const geometry::Atoms & atoms,
                           const basis::Basis & basis);

scf::FockBuilder<double>
generate_fock_builder(const basis::Basis & basis,
                      const arma::mat & one_electron_integral,
                      const arma::mat & two_electron_integral);


nlohmann::json rhf(const nlohmann::json & input,
                   const geometry::Atoms & atoms,
                   const basis::Basis & basis);
}

namespace hfincpp::hf::gradient {
arma::cube core_hamiltonian(const geometry::Atoms & atoms,
                            const basis::Basis & basis);

arma::cube two_electron_integral(const basis::Basis & basis);
}


#endif //HF_RHF_H

#ifndef HF_RHF_H
#define HF_RHF_H

#include <json.hpp>

#include "gradient/driver.h"
#include "geometry/geometry.h"
#include "basis/basis.h"
#include "scf/scf.h"
#include "scf/setup.h"

namespace hfincpp::hf {

template<typename T>
struct RHFSetup : scf::Setup<T> {
  arma::mat H0;
  scf::FockBuilder<T> fock_builder;
  scf::EnergyBuilder<T> energy_builder;
  scf::OccupationBuilder occupation_builder;
  scf::DensityMatrix<T> initial_guess;
};

template<typename T>
using RHFResult = scf::Result<RHFSetup<T>, T>;

arma::mat core_hamiltonian(const geometry::Atoms & atoms,
                           const basis::Basis & basis);

scf::FockBuilder<double>
generate_fock_builder(const basis::Basis & basis,
                      const arma::mat & one_electron_integral,
                      const arma::mat & two_electron_integral);

template<typename T>
nlohmann::json rhf(const RHFSetup<T> & setup);

nlohmann::json rhf(const nlohmann::json & input,
                   const geometry::Atoms & atoms,
                   const basis::Basis & basis);
}

namespace hfincpp::hf::gradient {
arma::cube core_hamiltonian(const geometry::Atoms & atoms,
                            const basis::Basis & basis);

arma::cube two_electron_integral(const basis::Basis & basis);

hfincpp::gradient::GradientDriver driver(const nlohmann::json & input,
                                         const geometry::Atoms & atoms,
                                         const basis::Basis & basis);

hfincpp::gradient::EnergyDriver energy_driver(const nlohmann::json & input,
                                              const geometry::Atoms & old_atoms,
                                              const basis::Basis & old_basis);
}


#endif //HF_RHF_H

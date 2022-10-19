#include "rhf.h"

#include "scf/scf.h"
#include "scf/occupation.h"

#include "integral/obara_saika.h"
#include "integral/rys_quadrature.h"

#include "mixing/simple_mixing.h"

namespace hfincpp {
namespace hf {

arma::mat core_hamiltonian(const geometry::Atoms & atoms,
                           const basis::Basis & basis) {
  const arma::mat kinetic = integral::obara_saika::kinetic_integral(basis);
  const arma::mat nuclear_attraction =
      integral::rys_quadrature::nuclear_attraction_integral(atoms, basis);

  return kinetic + nuclear_attraction;
}

scf::FockBuilder<double> generate_fock_builder(const geometry::Atoms & atoms,
                                               const basis::Basis & basis) {

  const arma::mat coulomb_integral =
      integral::rys_quadrature::electron_repulsive_integral(basis);

  const arma::mat overlap = integral::obara_saika::overlap_integral(basis);
  const arma::mat H0 = core_hamiltonian(atoms, basis);

  arma::mat exchange_integral(arma::size(coulomb_integral));

  const auto n_ao = basis.n_functions();

#pragma omp parallel for
  for(int i=0; i<n_ao; i++) {
    for(int j=0; j<n_ao; j++) {
      for(int k=0; k<n_ao; k++) {
        for(int l=0; l<n_ao; l++) {
          exchange_integral(i + l * n_ao, k + j * n_ao) =
              coulomb_integral(i + j * n_ao, k + l * n_ao);
        }
      }
    }
  }

  const arma::mat two_electron_integral = coulomb_integral - 0.5 * exchange_integral;

  return [n_ao, H0,
          two_electron_integral](const scf::DensityMatrix<double> & density) {
    scf::FockMatrix<double> result;
    for(const auto & i_density : density) {
      result.push_back(H0 +
      arma::reshape(two_electron_integral * arma::vectorise(i_density), n_ao, n_ao));
    }

    return result;
  };


}

nlohmann::json rhf(const nlohmann::json & input,
                   const geometry::Atoms & atoms,
                   const basis::Basis & basis) {

  const scf::OverlapMatrix<double> overlap = {integral::obara_saika::overlap_integral(basis)};
  const double n_elec_per_orb = 2.0;
  const scf::OccupationBuilder occupation_builder =
      scf::occupation::simple_occupation(n_elec_per_orb);

  const scf::FockBuilder<double> fock_builder = generate_fock_builder(atoms, basis);

  scf::DensityMatrix<double> initial_guess;

  const int n_elec = atoms.n_elec();

  const std::string initial_guess_method = input.at("initial_guess");
  if(initial_guess_method == "H0") {
    const arma::mat H0 = core_hamiltonian(atoms, basis);

    arma::cx_vec eigvals;
    arma::cx_mat eigvecs;
    arma::eig_pair(eigvals, eigvecs, H0, overlap[0]);

    assert(H0.is_hermitian() && overlap[0].is_hermitian());
    const arma::vec real_eigvals = arma::real(eigvals);
    const arma::uvec sort_index = arma::sort_index(real_eigvals);
    const arma::vec occupation_vector = occupation_builder(real_eigvals(sort_index), n_elec);
    const arma::mat sorted_orbitals = arma::real(eigvecs.cols(sort_index));
    initial_guess = {sorted_orbitals * arma::diagmat(occupation_vector) * sorted_orbitals.t()};
  } else {
    throw Error("only initial guess method of H0 is allowed for current version");
  }

  if(input.at("mixing") == "simple_mixing") {
    const double alpha = input.at("mixing_alpha");
    const auto update_method =
        mixing::simple_mixing<scf::DensityMatrix<double>>{alpha, initial_guess};

  }
}

}
}
#include "rhf.h"

#include "scf/scf.h"
#include "scf/occupation.h"

#include "integral/obara_saika.h"
#include "integral/rys_quadrature.h"

#include "mixing/simple_mixing.h"

namespace hfincpp::hf {

arma::mat core_hamiltonian(const geometry::Atoms & atoms,
                           const basis::Basis & basis) {
  const arma::mat kinetic = integral::obara_saika::kinetic_integral(basis);
  const arma::mat nuclear_attraction =
      integral::rys_quadrature::nuclear_attraction_integral(atoms, basis);

  return kinetic + nuclear_attraction;
}

scf::FockBuilder<double>
generate_fock_builder(const basis::Basis & basis,
                      const arma::mat & one_electron_integral,
                      const arma::mat & two_electron_integral) {


  const auto n_ao = basis.n_functions();
  const arma::mat overlap = integral::obara_saika::overlap_integral(basis);
  const arma::mat & H0 = one_electron_integral;

  return [n_ao, H0,
      two_electron_integral](const scf::DensityMatrix<double> & density) {
    scf::FockMatrix<double> result(arma::size(density));
    for (arma::uword i = 0; i < density.n_slices; i++) {
      result.slice(i) = H0 +
                        arma::reshape(
                            two_electron_integral *
                            arma::vectorise(density.slice(i)),
                            n_ao, n_ao);
    }

    return result;
  };


}

scf::EnergyBuilder<double>
generate_energy_builder(const geometry::Atoms & atoms,
                        const arma::mat & one_electron_integral,
                        const arma::mat & two_electron_integral) {

  const arma::mat & H0 = one_electron_integral;
  const arma::vec norm_squared = arma::sum(arma::square(atoms.xyz)).t();
  arma::mat distance_squared = -2.0 * atoms.xyz.t() * atoms.xyz;
  distance_squared.each_col() += norm_squared;
  distance_squared.each_row() += norm_squared.t();
  arma::mat inverse_r = 1.0 / arma::sqrt(distance_squared);
  inverse_r.diag().zeros();
  const double repulsion_among_cores =
      0.5 * arma::dot(atoms.atomic_numbers, inverse_r * atoms.atomic_numbers);

  return [H0, two_electron_integral, repulsion_among_cores](
      const scf::DensityMatrix<double> & density) {
    double result = 0;
    for (arma::uword i = 0; i < density.n_slices; i++) {
      result +=
          0.5 * arma::dot(arma::vectorise(density.slice(i)),
                          two_electron_integral *
                          arma::vectorise(density.slice(i)))
          + arma::dot(density.slice(i), H0);
    }

    return result + repulsion_among_cores;
  };

}

RHFSetup<double> setup(const nlohmann::json & input,
                       const geometry::Atoms & atoms,
                       const basis::Basis & basis) {

  const int print_level = input.at("print_level");
  double hf_timer_pool = 0;
  Timer hf_timer;

  scf::OverlapMatrix<double> overlap(basis.n_functions(), basis.n_functions(),
                                     1);
  overlap.slice(0) = integral::obara_saika::overlap_integral(basis);

  double overlap_build_time = hf_timer.elapsed() - hf_timer_pool;
  hf_timer_pool += overlap_build_time;

  if (print_level > 1) {
    fmt::print("overlap matrix build comsumed: {} s\n",
               format(overlap_build_time));
  }

  const arma::mat coulomb_integral =
      integral::rys_quadrature::electron_repulsive_integral(basis);

  arma::mat exchange_integral(arma::size(coulomb_integral));

  const auto n_ao = basis.n_functions();

#pragma omp parallel for
  for (int i = 0; i < n_ao; i++) {
    for (int j = 0; j < n_ao; j++) {
      for (int k = 0; k < n_ao; k++) {
        for (int l = 0; l < n_ao; l++) {
          exchange_integral(i + l * n_ao, k + j * n_ao) =
              coulomb_integral(i + j * n_ao, k + l * n_ao);
        }
      }
    }
  }

  const arma::mat two_electron_integral =
      coulomb_integral - 0.5 * exchange_integral;

  double two_electron_integral_build_time = hf_timer.elapsed() - hf_timer_pool;
  hf_timer_pool += two_electron_integral_build_time;
  if (print_level > 1) {
    fmt::print("2-electron integral build comsumed: {} s\n",
               format(two_electron_integral_build_time));
  }
  const arma::mat H0 = core_hamiltonian(atoms, basis);
  double one_electron_integral_build_time = hf_timer.elapsed() - hf_timer_pool;
  hf_timer_pool += one_electron_integral_build_time;
  if (print_level > 1) {
    fmt::print("1-electron integral build comsumed: {} s\n",
               format(one_electron_integral_build_time));
  }

  const double n_elec_per_orb = 2.0;
  const scf::OccupationBuilder occupation_builder =
      scf::occupation::simple_occupation(n_elec_per_orb);

  const scf::FockBuilder<double> fock_builder =
      generate_fock_builder(basis, H0, two_electron_integral);

  const scf::EnergyBuilder<double> energy_builder =
      generate_energy_builder(atoms, H0, two_electron_integral);

  RHFSetup<double> setup;
  setup.atoms = atoms;
  setup.overlap = overlap;
  setup.basis = basis;
  setup.occupation_builder = occupation_builder;
  setup.fock_builder = fock_builder;
  setup.energy_builder = energy_builder;
  setup.H0 = H0;

  return setup;
}

scf::DensityMatrix<double>
    initial_guess(const nlohmann::json & input,
                  const arma::mat & H0,
                  const scf::OverlapMatrix<double> overlap,
                  const scf::OccupationBuilder occupation_builder,
                  const int n_elec) {

  const std::string initial_guess_method = input.at("initial_guess");

  if (initial_guess_method == "H0") {
    scf::DensityMatrix<double> guess(arma::size(overlap));

    arma::cx_vec eigvals;
    arma::cx_mat eigvecs;

    arma::eig_pair(eigvals, eigvecs, H0, overlap.slice(0));

    assert(H0.is_hermitian() && overlap.slice(0).is_hermitian());
    const arma::vec real_eigvals = arma::real(eigvals);
    const arma::uvec sort_index = arma::sort_index(real_eigvals);
    const arma::vec occupation_vector = occupation_builder(
        real_eigvals(sort_index), n_elec);

    arma::mat sorted_orbitals = arma::real(eigvecs.cols(sort_index));

    const arma::rowvec normalization_constant =
        arma::real(
            arma::sum(sorted_orbitals % (overlap.slice(0) * sorted_orbitals)));

    sorted_orbitals.each_row() %= 1.0 / arma::sqrt(normalization_constant);
    const arma::mat density =
        sorted_orbitals * arma::diagmat(occupation_vector) *
        sorted_orbitals.t();

    guess.slice(0) = density;

    return guess;
  } else {
    throw Error(
        "only initial guess method of H0 is allowed for current version");
  }

}

RHFResult<double> rhf(const nlohmann::json & input,
                      const RHFSetup<double> & setup) {

  const int print_level = input.at("print_level");

  const int max_iter = input.at("max_iter");
  const double energy_tolerance = input.at("energy_tolerance");

  const scf::DensityMatrix<double>
  guess = initial_guess(input,
                        setup.H0,
                        setup.overlap,
                        setup.occupation_builder,
                        setup.atoms.n_elec());

  if (input.at("mixing") == "simple_mixing") {
    const double alpha = input.at("mixing_alpha");
    const auto update_method =
        mixing::simple_mixing<scf::DensityMatrix<double>>{alpha, guess};

    const auto scf_result = scf::scf<
        mixing::simple_mixing<scf::DensityMatrix<double>>,
        double>(setup.energy_builder,
                setup.fock_builder,
                setup.occupation_builder,
                update_method,
                setup.overlap,
                guess,
                {(double) setup.atoms.n_elec()},
                max_iter, energy_tolerance, print_level);

    return {setup, scf_result};
  } else {
    throw Error("Only simple mixing method is allowed in current method");
  }
}

nlohmann::json rhf(const nlohmann::json & input,
                   const geometry::Atoms & atoms,
                   const basis::Basis & basis) {
  const auto scf_setup = setup(input, atoms, basis);


  const auto result = rhf(input, scf_setup);
  auto scf_json = result.scf_result.to_json();
  util::put(scf_json, "basis_labels", basis.function_labels);
  return scf_json;
}

}

namespace hfincpp::hf::gradient {

arma::cube core_hamiltonian(const geometry::Atoms & atoms,
                            const basis::Basis & basis) {
  const arma::cube kinetic =
      integral::obara_saika::gradient::kinetic_integral(basis);
  const arma::cube nuclear_attraction =
      integral::rys_quadrature::gradient::nuclear_attraction_integral(atoms,
                                                                      basis);

  return kinetic + nuclear_attraction;
}

arma::cube two_electron_integral(const basis::Basis & basis) {

  const arma::cube coulomb_integral =
      integral::rys_quadrature::gradient::electron_repulsive_integral(basis);

  arma::cube exchange_integral(arma::size(coulomb_integral));

  const auto n_ao = basis.n_functions();

#pragma omp parallel for
  for (arma::uword slice = 0; slice < coulomb_integral.n_slices; slice++) {
    for (int i = 0; i < n_ao; i++) {
      for (int j = 0; j < n_ao; j++) {
        for (int k = 0; k < n_ao; k++) {
          for (int l = 0; l < n_ao; l++) {
            exchange_integral(i + l * n_ao, k + j * n_ao, slice) =
                coulomb_integral(i + j * n_ao, k + l * n_ao, slice);
          }
        }
      }
    }
  }


  return coulomb_integral - 0.5 * exchange_integral;

}

hfincpp::gradient::GradientDriver driver(const nlohmann::json & input,
                                         const geometry::Atoms & old_atoms,
                                         const basis::Basis & old_basis) {

  return
      [input, old_atoms, old_basis]
          (const geometry::Atoms & atoms) -> std::pair<double, arma::mat> {
        const basis::Basis basis(atoms, old_basis.basis_name);
        const auto scf_setup = setup(input, atoms, basis);
        const auto result = rhf(input, scf_setup);
        const arma::cube gradient_H0 = core_hamiltonian(atoms, basis);
        const arma::cube gradient_2e = two_electron_integral(basis);
        const arma::cube gradient_overlap =
            integral::obara_saika::gradient::overlap_integral(basis);

        const arma::mat density = result.scf_result.density.slice(0);
        const arma::vec vectorised_density = arma::vectorise(density);
        const arma::mat orbitals = result.scf_result.orbitals.slice(0);
        const arma::vec eigenvalues = result.scf_result.eigenvalues.col(0);
        const arma::vec occupations = result.scf_result.occupations.col(0);
        const arma::mat energy_weighted_density =
            orbitals * arma::diagmat(eigenvalues % occupations) * orbitals.t();

        arma::mat gradient(arma::size(atoms.xyz));
        for(arma::uword i=0; i<gradient.n_elem; i++) {
          gradient(i) =
              arma::accu(gradient_H0.slice(i) % density)
              - arma::accu(energy_weighted_density % gradient_overlap.slice(i))
              + 0.5 * arma::dot(vectorised_density,
                                gradient_2e.slice(i) * vectorised_density);

        }

        return {result.scf_result.energy, gradient};
      };
}

}
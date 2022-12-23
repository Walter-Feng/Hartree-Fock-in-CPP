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

nlohmann::json rhf(const nlohmann::json & input,
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

  scf::DensityMatrix<double> initial_guess(arma::size(overlap));

  const int n_elec = atoms.n_elec();

  const std::string initial_guess_method = input.at("initial_guess");
  const int max_iter = input.at("max_iter");
  const double energy_tolerance = input.at("energy_tolerance");

  if (initial_guess_method == "H0") {
    const arma::mat H0 = core_hamiltonian(atoms, basis);

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

    initial_guess.slice(0) = density;
  } else {
    throw Error(
        "only initial guess method of H0 is allowed for current version");
  }

  if (input.at("mixing") == "simple_mixing") {
    const double alpha = input.at("mixing_alpha");
    const auto update_method =
        mixing::simple_mixing<scf::DensityMatrix<double>>{alpha, initial_guess};

    const auto scf_result = scf::scf<
        mixing::simple_mixing<scf::DensityMatrix<double>>,
        double>(energy_builder, fock_builder, occupation_builder,
                update_method, overlap, initial_guess,
                {(double) atoms.n_elec()},
                max_iter, energy_tolerance, print_level);

    auto scf_json = scf_result.to_json();
    util::put(scf_json, "basis_labels", basis.function_labels);
    return scf_json;
  } else {
    throw Error("Only simple mixing method is allowed in current method");
  }
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

nlohmann::json rhf(const nlohmann::json & input,
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

  scf::DensityMatrix<double> initial_guess(arma::size(overlap));

  const int n_elec = atoms.n_elec();

  const std::string initial_guess_method = input.at("initial_guess");
  const int max_iter = input.at("max_iter");
  const double energy_tolerance = input.at("energy_tolerance");

  if (initial_guess_method == "H0") {
    const arma::mat H0 = core_hamiltonian(atoms, basis);

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

    initial_guess.slice(0) = density;
  } else {
    throw Error(
        "only initial guess method of H0 is allowed for current version");
  }

  if (input.at("mixing") == "simple_mixing") {
    const double alpha = input.at("mixing_alpha");
    const auto update_method =
        mixing::simple_mixing<scf::DensityMatrix<double>>{alpha, initial_guess};

    const auto scf_result = scf::scf<
        mixing::simple_mixing<scf::DensityMatrix<double>>,
        double>(energy_builder, fock_builder, occupation_builder,
                update_method, overlap, initial_guess,
                {(double) atoms.n_elec()},
                max_iter, energy_tolerance, print_level);

    auto scf_json = scf_result.to_json();
    util::put(scf_json, "basis_labels", basis.function_labels);
    return scf_json;
  } else {
    throw Error("Only simple mixing method is allowed in current method");
  }
}

}
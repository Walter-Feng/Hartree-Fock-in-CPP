#include "catch.h"

#include "rhf.h"

#include "gradient/numerical.h"
#include "integral/rys_quadrature.h"

TEST_CASE("Check RHF Implementation") {

  using namespace hfincpp;

  SECTION("Check core hamiltonian gradient") {

    hfincpp::geometry::Atoms atoms;
    atoms.atomic_numbers = {1, 9};
    atoms.xyz = {
        {0, 3.77945},
        {0, 0},
        {0, 0}
    };

    atoms.symbols = {"H", "F"};

    const std::string basis_name = "6-31g";
    hfincpp::basis::Basis basis(atoms, basis_name);
    const std::function<arma::mat(const hfincpp::geometry::Atoms &)>
        numerical_gradient_functor =
        [basis_name](const hfincpp::geometry::Atoms & atoms) -> arma::mat {
          const hfincpp::basis::Basis basis(atoms, basis_name);

          return hf::core_hamiltonian(atoms, basis);
        };

    const arma::cube analytical = hf::gradient::core_hamiltonian(atoms, basis);

    const auto numerical_gradient =
        hfincpp::gradient::numerical(numerical_gradient_functor, atoms);

    arma::cube numerical_in_cube(arma::size(analytical));
    for (arma::uword i = 0; i < numerical_in_cube.n_slices; i++) {
      numerical_in_cube.slice(i) = numerical_gradient[i];
    }

    CHECK(arma::abs(analytical - numerical_in_cube).max() < 1e-8);

  }


  SECTION("Check fock matrix gradient") {

    hfincpp::geometry::Atoms atoms;
    atoms.atomic_numbers = {1, 9};
    atoms.xyz = {
        {0, 3.77945},
        {0, 0},
        {0, 0}
    };

    atoms.symbols = {"H", "F"};

    const std::string basis_name = "6-31g";
    hfincpp::basis::Basis basis(atoms, basis_name);
    const scf::DensityMatrix<double> density =
        arma::ones(basis.n_functions(), basis.n_functions(), 1);
    const std::function<arma::mat(const hfincpp::geometry::Atoms &)>
        numerical_gradient_functor =
        [basis_name, density](
            const hfincpp::geometry::Atoms & atoms) -> arma::mat {
          const hfincpp::basis::Basis basis(atoms, basis_name);
          const arma::mat H0 = hf::core_hamiltonian(atoms, basis);

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

          const auto fock_builder =
              hf::generate_fock_builder(basis, H0, two_electron_integral);

          return fock_builder(density);
        };

    const auto numerical_gradient =
        hfincpp::gradient::numerical(numerical_gradient_functor, atoms);

    arma::cube analytical = hf::gradient::core_hamiltonian(atoms, basis);
    const arma::cube two_electron_integral =
        hf::gradient::two_electron_integral(basis);

    for (arma::uword i = 0; i < analytical.n_slices; i++) {
      analytical.slice(i) +=
          arma::reshape(two_electron_integral.slice(i) *
                        arma::vectorise(density.slice(0)),
                        basis.n_functions(),
                        basis.n_functions());
    }

    arma::cube numerical_in_cube(arma::size(analytical));
    for (arma::uword i = 0; i < numerical_in_cube.n_slices; i++) {
      numerical_in_cube.slice(i) = numerical_gradient[i];
    }

    CHECK(arma::abs(analytical - numerical_in_cube).max() < 1e-8);

  }
}
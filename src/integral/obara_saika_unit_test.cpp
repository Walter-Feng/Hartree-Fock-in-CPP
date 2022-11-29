#include "catch.h"

#include "obara_saika.h"

using namespace hfincpp::integral::obara_saika;

TEST_CASE("Check Basis struct") {
  SECTION("check file parse") {
    hfincpp::geometry::Atoms atoms;
    atoms.atomic_numbers = {1, 1};
    atoms.xyz = {
        {0, 3.77945},
        {0, 0},
        {0, 0}
    };

    atoms.symbols = {"H", "H"};

    hfincpp::basis::Basis basis(atoms, "cc-pvdz");

    const arma::mat overlap = overlap_integral(basis);
    const arma::mat kinetic = kinetic_integral(basis);
    const arma::mat eri = electron_repulsive_integral(basis);

  }


  SECTION("Electron Repulsive Integral, ERI") {
    hfincpp::integral::GaussianFunction A;
    A.center = {0.0, 0.0, 0.0};
    A.angular = {0, 1, 1};
    A.exponent = 1.0;
    A.coef = 1.0;

    hfincpp::integral::GaussianFunction B;
    B.center = {0.5, 0.6, 0.7};
    B.angular = {1, 2, 0};
    B.exponent = 2.0;
    B.coef = 1.0;

    hfincpp::integral::GaussianFunction C;
    C.center = {0.9, 1.0, 1.1};
    C.angular = {2, 1, 1};
    C.exponent = 3.0;
    C.coef = 1.0;

    hfincpp::integral::GaussianFunction D;
    D.center = {1.0, 0.9, 0.8};
    D.angular = {1, 1, 0};
    D.exponent = 4.0;
    D.coef = 1.0;

    hfincpp::integral::ERI eri{A, B, C, D};

    const double eri_integral = electron_repulsive_integral(eri);

    CHECK(std::abs(eri_integral + 7.486283316104355e-08) < 1e-12);
  }

}
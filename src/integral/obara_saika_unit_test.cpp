#include <catch.hpp>

#include "obara_saika.h"

using namespace hfincpp::integral::obara_saika;

TEST_CASE("Check Basis struct") {
  SECTION("check file parse") {
    hfincpp::geometry::Atoms atoms;
    atoms.atomic_numbers = {8, 1, 1};
    atoms.xyz = {
        {0, 0, 0},
        {0, 1, 0},
        {0, 0, -1}
    };

    atoms.symbols = {"O", "H", "H"};

    hfincpp::basis::Basis basis(atoms, "cc-pvdz");

    const arma::mat overlap = overlap_integral(basis);
    const arma::mat kinetic = kinetic_integral(basis);

  }
}
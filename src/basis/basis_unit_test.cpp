#include "catch.h"

#include "basis.h"

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

    hfincpp::basis::Basis basis(atoms, "3-21g");
  }
}
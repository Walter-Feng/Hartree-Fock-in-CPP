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

  SECTION("arma playground") {
    arma::Mat<int> field = arma::zeros<arma::Mat<int>>(5,5);
    for(int i=0; i<5; i++) {
      for(int j=i; j<5; j++) {
        field(i, j) = (5 + 5 - i + 1) * i / 2 + j - i;
      }
    }

    field.print();
  }
}
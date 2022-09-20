#include <catch.hpp>

#include "rys_quadrature.h"

TEST_CASE("Check Rys quadrature integral implementation") {

  SECTION("Horizontal Recursion Relation, HRR") {
    using namespace integral::rys_quadrature;
    IntegralInfo info;
    info.polynomial = RysPolynomial{arma::vec{1.0}};
    info.a = 3;
    info.b = 4;
    info.c = 2;
    info.d = 5;

    info.p = 1;
    info.P = 0;
    info.q = 1;
    info.Q = 1;

    info.A = 0.0;
    info.B = 1.5;
    info.C = 2.0;
    info.D = 5.0;

    const auto result = horizontal_recursion_relation({info});

    CHECK(result.size() == std::pow(2, info.b) * std::pow(2, info.d));
    for(auto info_i : result) {
      CHECK(info_i.b == 0);
      CHECK(info_i.d == 0);

      CHECK(info_i.polynomial.coef.n_elem == 1);
    }
  }

  SECTION("Vertical Recursion Relation, VRR") {
    using namespace integral::rys_quadrature;
    IntegralInfo info;
    info.polynomial = RysPolynomial{arma::vec{1.0}};
    info.a = 3;
    info.b = 0;
    info.c = 2;
    info.d = 0;

    info.p = 1;
    info.P = 6;
    info.q = 1;
    info.Q = 9;

    info.A = 0.0;
    info.B = 1.5;
    info.C = 2.0;
    info.D = 5.0;

    const auto result = vertical_recursion_relation({info});

    for(const auto & i_info : result) {
      CHECK(i_info.a == 0);
      CHECK(i_info.b == 0);
      CHECK(i_info.c == 0);
      CHECK(i_info.d == 0);
    }
  }

  SECTION("Reduce to Rys Polynomial") {
    using namespace integral::rys_quadrature;
    IntegralInfo info;
    info.polynomial = RysPolynomial{arma::vec{1.0}};
    info.a = 3;
    info.b = 4;
    info.c = 2;
    info.d = 5;

    info.p = 1;
    info.P = 6;
    info.q = 1;
    info.Q = 9;

    info.A = 0.0;
    info.B = 1.5;
    info.C = 2.0;
    info.D = 5.0;

    const auto result = reduce_to_rys_polynomial(info);

    std::cout << result.coef << std::endl;
  }


}
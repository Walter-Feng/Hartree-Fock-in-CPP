#include <catch.hpp>

#include "rys_quadrature.h"
#include "obara_saika.h"

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
    for (auto info_i: result) {
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

    for (const auto & i_info: result) {
      CHECK(i_info.a == 0);
      CHECK(i_info.b == 0);
      CHECK(i_info.c == 0);
      CHECK(i_info.d == 0);
    }
  }

//  SECTION("Reduce to Rys Polynomial") {
//    using namespace integral::rys_quadrature;
//    IntegralInfo info;
//    info.polynomial = RysPolynomial{arma::vec{1.0}};
//    info.a = 3;
//    info.b = 4;
//    info.c = 2;
//    info.d = 5;
//
//    info.p = 1;
//    info.P = 6;
//    info.q = 1;
//    info.Q = 9;
//
//    info.A = 0.0;
//    info.B = 1.5;
//    info.C = 2.0;
//    info.D = 5.0;
//
//    const auto result = reduce_to_rys_polynomial(info);
//  }


  SECTION("Reduce to Rys Polynomial") {
    using namespace integral::rys_quadrature;

    integral::GaussianFunction A;
    A.center = {0.0, 0.0, 0.0};
    A.angular = {0, 1, 1};
    A.exponent = 1.0;
    A.coef = 1.0;

    integral::GaussianFunction B;
    B.center = {0.5, 0.6, 0.7};
    B.angular = {1, 2, 0};
    B.exponent = 2.0;
    B.coef = 1.0;

    integral::GaussianFunction C;
    C.center = {0.9, 1.0, 1.1};
    C.angular = {2, 1, 1};
    C.exponent = 3.0;
    C.coef = 1.0;

    integral::GaussianFunction D;
    D.center = {1.0, 0.9, 0.8};
    D.angular = {1, 1, 0};
    D.exponent = 4.0;
    D.coef = 1.0;

    integral::ERI eri{A, B, C, D};

    integral::ERI eri_prime{A, B, C, D};

    const double eri_integral = electron_repulsive_integral(eri);

    CHECK(std::abs(eri_integral + 7.486283316104355e-08) < 1e-12);
  }
}
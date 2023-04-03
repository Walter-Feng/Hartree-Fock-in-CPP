#include "catch.h"

#include "rys_quadrature.h"
#include "../rys_quadrature.h"

TEST_CASE("ERI kernel") {
  using namespace hfincpp::integral::faster_kernel::rys_quadrature;

  SECTION("electron_repulsive_integral_kernel") {

    double C00 = arma::randu();
    double D00 = arma::randu();
    double B10 = arma::randu();
    double B01 = arma::randu();
    double B00 = arma::randu();
    arma::vec buffer(20);
    int a = 3;
    int b = 3;

    electron_repulsive_integral_kernel<3, 3, 1>(buffer.memptr() + 4,
                                                C00, D00, B01, B10, B00);

    double analytical_value =
        6 * std::pow(B00, 3) + 18 * std::pow(B00, 2) * C00 * D00 +
        9 * B00 * (B10 +
                   std::pow(C00, 2)) * (B01 + std::pow(D00, 2)) +
        C00 * (3 * B10 +
               std::pow(C00, 2)) * D00 * (3 * B01 + std::pow(D00, 2));

    CHECK(buffer((a + 1) * (b + 2) - 1) == analytical_value);

    CHECK(binomial_coefficient(4, 2) == 6);
    CHECK(binomial_coefficient(6, 2) == 15);

  }


  SECTION("electron_repulsive_integral_conventional_kernel") {
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

    const double eri_integral =
        to_conventional::electron_repulsive_integral<5, 6>(eri);

    CHECK(std::abs(eri_integral + 7.486283316104355e-08) < 1e-12);
  }

}
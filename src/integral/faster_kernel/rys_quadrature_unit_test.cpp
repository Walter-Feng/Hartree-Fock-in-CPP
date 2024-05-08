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
    arma::mat buffer(4,4);
    int a = 3;
    int b = 3;

    two_center_recursive_kernel<double, 3, 3>(buffer.memptr(), C00, D00, B01,
                                              B10, B00, 1);

    double analytical_value =
        6 * std::pow(B00, 3) + 18 * std::pow(B00, 2) * C00 * D00 +
        9 * B00 * (B10 +
                   std::pow(C00, 2)) * (B01 + std::pow(D00, 2)) +
        C00 * (3 * B10 +
               std::pow(C00, 2)) * D00 * (3 * B01 + std::pow(D00, 2));

    CHECK(std::abs(buffer((a + 1) * (b + 1) - 1) - analytical_value) < 1e-12);

    CHECK(binomial_coefficient(4, 2) == 6);
    CHECK(binomial_coefficient(6, 2) == 15);

  }


  SECTION("electron_repulsive_integral_conventional_kernel") {
    hfincpp::integral::GaussianFunction A;
    A.center = {0.0, 0.0, 0.0};
    A.angular = {1, 0, 0};
    A.exponent = 1.0;
    A.coef = 1.0;

    hfincpp::integral::GaussianFunction B;
    B.center = {0.5, 0.6, 0.7};
    B.angular = {1, 0, 0};
    B.exponent = 2.0;
    B.coef = 1.0;

    hfincpp::integral::GaussianFunction C;
    C.center = {0.9, 1.0, 1.1};
    C.angular = {0,1,0};
    C.exponent = 3.0;
    C.coef = 1.0;

    hfincpp::integral::GaussianFunction D;
    D.center = {1.0, 0.9, 0.8};
    D.angular = {0, 0, 1};
    D.exponent = 4.0;
    D.coef = 1.0;

    hfincpp::integral::ERI eri{A, B, C, D};

    const double eri_integral =
        to_conventional::electron_repulsive_integral<2, 2>(eri);

    const double eri_intergral_from_four_center_kernel =
        to_conventional::electron_repulsive_integral_from_four_center_kernel<1, 1, 1, 1>(
            eri);

    CHECK(std::abs(eri_intergral_from_four_center_kernel - eri_integral) < 1e-12);
  }

}
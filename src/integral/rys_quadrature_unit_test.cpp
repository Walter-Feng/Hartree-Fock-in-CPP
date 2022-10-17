#include <catch.hpp>

#include "rys_quadrature.h"

using namespace hfincpp::integral::rys_quadrature;

const auto numerical_gradient_functor = [](
    const std::function<double(double)> & functor,
    const double original,
    const double step = 1e-3) -> double {
  return (functor(original + step) - functor(original - step)) / 2.0 / step;
};

double electron_repulsive_integral_numerical_gradient(
    const hfincpp::integral::ERI & info,
    const arma::Mat<int>::fixed<3, 4> & derivative_operator
) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      if (derivative_operator(i, j) > 0) {
        arma::Mat<int>::fixed<3, 4> new_operator = derivative_operator;
        new_operator(i, j)--;
        hfincpp::integral::ERI info_copy = info;
        const auto index = i + j * 3;
        const auto eri_functor = [&info_copy, index, &new_operator](
            double coord_value) -> double {

          switch (index) {
            case 0:
              info_copy.A.center(0) = coord_value;
              break;
            case 1:
              info_copy.A.center(1) = coord_value;
              break;
            case 2:
              info_copy.A.center(2) = coord_value;
              break;
            case 3:
              info_copy.B.center(0) = coord_value;
              break;
            case 4:
              info_copy.B.center(1) = coord_value;
              break;
            case 5:
              info_copy.B.center(2) = coord_value;
              break;
            case 6:
              info_copy.C.center(0) = coord_value;
              break;
            case 7:
              info_copy.C.center(1) = coord_value;
              break;
            case 8:
              info_copy.C.center(2) = coord_value;
              break;
            case 9:
              info_copy.D.center(0) = coord_value;
              break;
            case 10:
              info_copy.D.center(1) = coord_value;
              break;
            case 11:
              info_copy.D.center(2) = coord_value;
              break;
            default:
              break;
          }

          return electron_repulsive_integral_numerical_gradient(info_copy,
                                                                new_operator);
        };

        double original;
        switch (index) {
          case 0:
            original = info_copy.A.center(0);
            break;
          case 1:
            original = info_copy.A.center(1);
            break;
          case 2:
            original = info_copy.A.center(2);
            break;
          case 3:
            original = info_copy.B.center(0);
            break;
          case 4:
            original = info_copy.B.center(1);
            break;
          case 5:
            original = info_copy.B.center(2);
            break;
          case 6:
            original = info_copy.C.center(0);
            break;
          case 7:
            original = info_copy.C.center(1);
            break;
          case 8:
            original = info_copy.C.center(2);
            break;
          case 9:
            original = info_copy.D.center(0);
            break;
          case 10:
            original = info_copy.D.center(1);
            break;
          case 11:
            original = info_copy.D.center(2);
            break;
          default:
            original = 0;
            break;
        }

        return numerical_gradient_functor(eri_functor, original);
      }
    }
  }

  return electron_repulsive_integral(info);
}

TEST_CASE("Check Rys quadrature integral implementation") {

  SECTION("Horizontal Recursion Relation, HRR") {

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

  SECTION("Gradient of ERI") {
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

    arma::Mat<int>::fixed<3, 4> derivative_operator = {
        {1, 0, 0, 1},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };

    double diff =
        gradient::electron_repulsive_integral(eri, derivative_operator) -
        electron_repulsive_integral_numerical_gradient(eri, derivative_operator);

    CHECK(std::abs(diff) < 1e-10);
  }
}
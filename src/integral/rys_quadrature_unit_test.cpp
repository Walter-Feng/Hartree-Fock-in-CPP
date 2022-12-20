#include "catch.h"

#include "rys_quadrature.h"

#include "gradient/numerical.h"

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

        double original = 0;
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
            break;
        }

        return numerical_gradient_functor(eri_functor, original);
      }
    }
  }

  return electron_repulsive_integral(info);
}


double nuclear_attraction_integral_numerical_gradient(
    const hfincpp::integral::GaussianFunctionPair & pair,
    const arma::vec3 & center,
    const double charge,
    const arma::Mat<int>::fixed<3, 3> & derivative_operator
) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (derivative_operator(i, j) > 0) {
        arma::Mat<int>::fixed<3, 3> new_operator = derivative_operator;
        new_operator(i, j)--;
        auto pair_copy = pair;
        auto center_copy = center;
        const auto index = i + j * 3;
        const auto eri_functor =
            [&pair_copy,
             &center_copy,
             index,
             &new_operator,
             charge](double coord_value) -> double {

          switch (index) {
            case 0:
              pair_copy.first.center(0) = coord_value;
              break;
            case 1:
              pair_copy.first.center(1) = coord_value;
              break;
            case 2:
              pair_copy.first.center(2) = coord_value;
              break;
            case 3:
              pair_copy.second.center(0) = coord_value;
              break;
            case 4:
              pair_copy.second.center(1) = coord_value;
              break;
            case 5:
              pair_copy.second.center(2) = coord_value;
              break;
            case 6:
              center_copy(0) = coord_value;
              break;
            case 7:
              center_copy(1) = coord_value;
              break;
            case 8:
              center_copy(2) = coord_value;
              break;
            default:
              break;
          }

          return nuclear_attraction_integral_numerical_gradient(pair_copy,
                                                                center_copy,
                                                                charge,
                                                                new_operator);
        };

        double original = 0;
        switch (index) {
          case 0:
            original = pair_copy.first.center(0);
            break;
          case 1:
            original = pair_copy.first.center(1);
            break;
          case 2:
            original = pair_copy.first.center(2);
            break;
          case 3:
            original = pair_copy.second.center(0);
            break;
          case 4:
            original = pair_copy.second.center(1);
            break;
          case 5:
            original =pair_copy.second.center(2);
            break;
          case 6:
            original = center_copy(0);
            break;
          case 7:
            original = center_copy(1);
            break;
          case 8:
            original = center_copy(2);
            break;
          default:
            break;
        }

        return numerical_gradient_functor(eri_functor, original);
      }
    }
  }

  return nuclear_attraction_integral(pair, center, charge);
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

  SECTION("Gradient of NAI") {
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

    const arma::vec3 core_center = {0.25, 0.25, 0.25};

    arma::Mat<int>::fixed<3, 3> derivative_operator = {
        {1, 0, 0},
        {1, 0, 0},
        {1, 0, 0}
    };

    std::cout << "analytical: " << gradient::nuclear_attraction_integral({A,B}, core_center, 1, derivative_operator)
    << " numerical: " << nuclear_attraction_integral_numerical_gradient({A,B}, core_center, 1, derivative_operator) << std::endl;
    double diff =
        gradient::nuclear_attraction_integral({A,B}, core_center, 1, derivative_operator) -
        nuclear_attraction_integral_numerical_gradient({A,B}, core_center, 1, derivative_operator);

    CHECK(std::abs(diff) < 1e-9);
  }

//
//  SECTION("Check gradient of ERI in actual systems") {
//    hfincpp::geometry::Atoms atoms;
//    atoms.atomic_numbers = {1, 9};
//    atoms.xyz = {
//        {0, 3.77945},
//        {0, 0},
//        {0, 0}
//    };
//
//    atoms.symbols = {"H", "F"};
//
//    const std::string basis_name = "6-31g";
//    hfincpp::basis::Basis basis(atoms, basis_name);
//
//    arma::Mat<int>::fixed<3, 4> gradient_operator_on_i = {
//        {1, 0, 0, 0},
//        {0, 0, 0, 0},
//        {0, 0, 0, 0}
//    };
//
//    arma::Mat<int>::fixed<3, 4> gradient_operator_on_j = {
//        {0, 1, 0, 0},
//        {0, 0, 0, 0},
//        {0, 0, 0, 0}
//    };
//
//    const arma::mat gradient_on_i =
//        gradient::electron_repulsive_integral(basis, gradient_operator_on_i);
//
//    const arma::mat gradient_on_j =
//        gradient::electron_repulsive_integral(basis, gradient_operator_on_j);
//
//    CHECK(arma::abs(gradient_on_j
//    - gradient::transpose_electron_repulsive_integral_i_with_j(gradient_on_i))
//      .max() < 1e-8);
//
//    const arma::cube gradient_atomic =
//        gradient::electron_repulsive_integral(basis);
//
//    const std::function<arma::mat(const hfincpp::geometry::Atoms &)>
//        numerical_eri_functor =
//        [basis_name](const hfincpp::geometry::Atoms & atoms) -> arma::mat {
//          const hfincpp::basis::Basis basis(atoms, basis_name);
//
//          return electron_repulsive_integral(basis);
//        };
//
//    const auto numerical_eri_gradient =
//        hfincpp::gradient::numerical(numerical_eri_functor, atoms);
//
//    arma::cube numerical_in_cube(arma::size(gradient_atomic));
//    for(arma::uword i=0; i<numerical_in_cube.n_slices; i++) {
//      numerical_in_cube.slice(i) = numerical_eri_gradient[i];
//    }
//
//    CHECK(arma::abs(gradient_atomic - numerical_in_cube).max() < 1e-8);
//
//  }

}
#include "catch.h"

#include "rys_quadrature.h"

TEST_CASE("ERI kernel") {
  using namespace hfincpp::integral::faster_kernel::rys_quadrature;
  double C00 = arma::randu();
  double D00 = arma::randu();
  double B10 = arma::randu();
  double B01 = arma::randu();
  double B00 = arma::randu();
  arma::vec buffer(20);
  int a = 3;
  int b = 3;

  electron_repulsive_integral_kernel<3,3,1>(buffer.memptr() + 4,
                                            C00,D00,B01,B10,B00);

  double analytical_value =
      6*std::pow(B00,3) + 18*std::pow(B00,2)*C00*D00 + 9*B00*(B10 +
        std::pow(C00,2))*(B01 + std::pow(D00,2)) + C00*(3*B10 +
        std::pow(C00,2))*D00*(3*B01 + std::pow(D00,2));

  CHECK(buffer((a+1) * (b+2) - 1) == analytical_value);

}
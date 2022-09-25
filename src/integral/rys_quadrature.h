#ifndef INTEGRAL_RYS_QUADRATURE_H
#define INTEGRAL_RYS_QUADRATURE_H

#include "integral.h"

namespace integral {
namespace rys_quadrature {

struct RysPolynomial {
  arma::vec coef;

  double operator()(double t_square) const;
  RysPolynomial operator*(double factor) const;
  RysPolynomial operator*(const RysPolynomial & another) const;

};

struct IntegralInfo {
  RysPolynomial polynomial;

  int a;
  int b;
  int c;
  int d;


  double p;
  double P;

  double q;
  double Q;

  double A;
  double B;
  double C;
  double D;

  double alpha; // exponent of A
  double beta; // exponent of B
  double gamma; // exponent of C
  double delta; // exponent of D

  int grad_a = 0;
  int grad_b = 0;
  int grad_c = 0;
  int grad_d = 0;

  [[nodiscard]] std::vector<IntegralInfo>
  horizontal_recursion_relation_b_to_a() const;

  [[nodiscard]] std::vector<IntegralInfo>
  horizontal_recursion_relation_c_to_d() const;

  [[nodiscard]] std::vector<IntegralInfo>
  vertical_recursion_relation_a() const;

  [[nodiscard]] std::vector<IntegralInfo>
  vertical_recursion_relation_c() const;
};


std::vector<IntegralInfo> horizontal_recursion_relation(const std::vector<IntegralInfo> & info);
std::vector<IntegralInfo> vertical_recursion_relation(const std::vector<IntegralInfo> & info);
double electron_repulsive_integral(const ERI & eri_info);

namespace gradient {
double electron_repulsive_integral(const ERI & eri_info,
                                   const arma::Mat<int>::fixed<3, 4> & derivative_operator);
}
}
}

#endif //INTEGRAL_RYS_QUADRATURE_H

#ifndef INTEGRAL_RYS_QUADRATURE_H
#define INTEGRAL_RYS_QUADRATURE_H

#include <armadillo>

namespace integral {
namespace rys_quadrature {

struct RysPolynomial {
  arma::vec coef;

  double operator()(double t_square) const;
  RysPolynomial operator*(double factor) const;
  RysPolynomial operator*(const RysPolynomial & another) const;

};


struct ERI {
  arma::vec3 A_coord;
  arma::vec3 B_coord;
  arma::vec3 C_coord;
  arma::vec3 D_coord;

  arma::Col<int> A_angular;
  arma::Col<int> B_angular;
  arma::Col<int> C_angular;
  arma::Col<int> D_angular;

  double A_exponent;
  double B_exponent;
  double C_exponent;
  double D_exponent;
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
RysPolynomial reduce_to_rys_polynomial(const IntegralInfo & info);
double electron_repulsive_integral(const ERI & eri_info);
}
}

#endif //INTEGRAL_RYS_QUADRATURE_H

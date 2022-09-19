#include "rys_quadrature.h"

#include <armadillo>

namespace integral {
namespace rys_quadrature {

struct RysPolynomial {
  arma::vec coef;

  double operator()(const double t_square) const {
    arma::vec powered(coef.n_elem);
    for (int i = 0; i < coef.n_elem; i++) {
      powered(i) = std::pow(t_square, i);
    }

    return arma::dot(coef, powered);
  }

  RysPolynomial operator*(const double factor) const {
    return {coef * factor};
  }

  RysPolynomial operator*(const RysPolynomial & another) const {
    arma::vec new_coef(another.coef.n_elem + this->coef.n_elem - 1,
                       arma::fill::zeros);

    for (int i = 0; i < another.coef.n_elem; i++) {
      new_coef(arma::span(i, coef.n_elem + i - 1)) += coef * another(i);
    }

    return {new_coef};
  }

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
  horizontal_recursion_relation_b_to_a() const {
    if (b == 0) {
      return {};
    } else {
      auto term_1 = *this;
      auto term_2 = *this;
      term_1.a++;
      term_1.b--;
      term_2.polynomial = term_2.polynomial * (term_2.A - term_2.B);
      term_2.b--;

      return {term_1, term_2};
    }
  }

  [[nodiscard]] std::vector<IntegralInfo>
  horizontal_recursion_relation_c_to_d() const {
    if (d == 0) {
      return {};
    } else {
      auto term_1 = *this;
      auto term_2 = *this;
      term_1.c ++;
      term_1.d --;
      term_2.polynomial = term_2.polynomial * (term_2.C - term_2.D);
      term_2.d --;

      return {term_1, term_2};
    }
  }
};

std::vector<IntegralInfo> horizontal_recursion_relation(const std::vector<IntegralInfo> & info) {

  std::vector<IntegralInfo> result = info;

  for(size_t i=0; i<info.size(); i++) {
    const auto & i_info = info[i];
    if(i_info.b > 0) {
      const auto b_horiz = i_info.horizontal_recursion_relation_b_to_a();

      const auto b_recursive = horizontal_recursion_relation(b_horiz);
    }
  }

}
}
}

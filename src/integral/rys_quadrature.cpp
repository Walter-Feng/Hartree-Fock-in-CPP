#include "rys_quadrature.h"
#include <cassert>
#include <armadillo>

extern "C" {
#include <rys_roots.h>
}

namespace integral::rys_quadrature {

double RysPolynomial::operator()(const double t_square) const {
  arma::vec powered(coef.n_elem);
  for (arma::uword i = 0; i < coef.n_elem; i++) {
    powered(i) = std::pow(t_square, i);
  }

  return arma::dot(coef, powered);
}

RysPolynomial RysPolynomial::operator*(const double factor) const {
  return {coef * factor};
}

RysPolynomial RysPolynomial::operator*(const RysPolynomial & another) const {
  arma::vec new_coef(another.coef.n_elem + this->coef.n_elem - 1,
                     arma::fill::zeros);

  for (arma::uword i = 0; i < another.coef.n_elem; i++) {
    new_coef(arma::span(i, coef.n_elem + i - 1)) += coef * another.coef(i);
  }

  return {new_coef};
}


std::vector<IntegralInfo>
IntegralInfo::horizontal_recursion_relation_b_to_a() const {
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
IntegralInfo::horizontal_recursion_relation_c_to_d() const {
  if (d == 0) {
    return {};
  } else {
    auto term_1 = *this;
    auto term_2 = *this;
    term_1.c++;
    term_1.d--;
    term_2.polynomial = term_2.polynomial * (term_2.C - term_2.D);
    term_2.d--;

    return {term_1, term_2};
  }
}

std::vector<IntegralInfo>
IntegralInfo::vertical_recursion_relation_a() const {
  if (a <= 0) {
    return {};
  } else if (a == 1) {
    // The second term, a * B_10 * I(a-1, 0, c, 0, t) is eliminated
    auto term_1 = *this;
    const RysPolynomial C_00 = {arma::vec{term_1.P - term_1.A,
                                          -term_1.q * (term_1.P - term_1.Q) /
                                          (term_1.p + term_1.q)}};

    term_1.a--;
    term_1.polynomial = term_1.polynomial * C_00;

    if (c <= 0) {
      // The third term, c * B_00 * I(a, 0, c-1, 0, t) is also eliminated
      return {term_1};
    } else {
      auto term_3 = *this;
      const RysPolynomial B_00 = {arma::vec{0, 0.5 / (term_3.p + term_3.q)}};
      term_3.a--;
      term_3.polynomial = term_3.polynomial * B_00 * term_3.c;
      term_3.c--;

      return {term_1, term_3};
    }

  } else {
    auto term_1 = *this;
    const RysPolynomial C_00 = {arma::vec{term_1.P - term_1.A,
                                          -term_1.q * (term_1.P - term_1.Q) /
                                          (term_1.p + term_1.q)}};

    term_1.a--;
    term_1.polynomial = term_1.polynomial * C_00;

    auto term_2 = *this;
    const RysPolynomial B_10 = {
        arma::vec{0.5 / term_2.p,
                  -0.5 * term_2.q / term_2.p / (term_2.p + term_2.q)}};

    term_2.polynomial = term_2.polynomial * B_10 * (term_2.a - 1);
    term_2.a -= 2;

    if (c == 0) {
      // The third term, c * B_00 * I(a, 0, c-1, 0, t) is eliminated
      return {term_1, term_2};
    } else {
      auto term_3 = *this;
      const RysPolynomial B_00 = {arma::vec{0, 0.5 / (term_3.p + term_3.q)}};

      term_3.a--;
      term_3.polynomial = term_3.polynomial * B_00 * term_3.c;
      term_3.c--;

      return {term_1, term_2, term_3};
    }
  }
}

[[nodiscard]] std::vector<IntegralInfo>
IntegralInfo::vertical_recursion_relation_c() const {
  if (c <= 0) {
    return {};
  } else if (c == 1) {
    // The third term, c * B_01 * I (a, 0, c-1, 0, t) is eliminated
    auto term_1 = *this;
    const RysPolynomial D_00 = {arma::vec{term_1.Q - term_1.C,
                                          +term_1.p * (term_1.P - term_1.Q) /
                                          (term_1.p + term_1.q)}};

    term_1.c--;
    term_1.polynomial = term_1.polynomial * D_00;

    if (a <= 0) {
      return {term_1};
    } else {
      auto term_2 = *this;
      const RysPolynomial B_00 = {arma::vec{0, 0.5 / (term_1.p + term_1.q)}};

      term_2.c--;
      term_2.polynomial = term_2.polynomial * B_00 * term_2.a;
      term_2.a--;

      return {term_1, term_2};
    }
  } else {
    // The third term, c * B_01 * I (a, 0, c-1, 0, t) is eliminated
    auto term_1 = *this;
    const RysPolynomial D_00 = {arma::vec{term_1.Q - term_1.C,
                                          +term_1.p * (term_1.P - term_1.Q) /
                                          (term_1.p + term_1.q)}};

    term_1.c--;
    term_1.polynomial = term_1.polynomial * D_00;

    auto term_3 = *this;
    const RysPolynomial B_01 = {arma::vec{0.5 / term_3.q,
                                          -0.5 * term_3.p / term_3.q /
                                          (term_3.p + term_3.q)}};

    term_3.polynomial = term_3.polynomial * B_01 * (term_3.c - 1);
    term_3.c -= 2;

    if (a <= 0) {
      return {term_1, term_3};
    } else {
      auto term_2 = *this;
      const RysPolynomial B_00 = {arma::vec{0, 0.5 / (term_2.p + term_2.q)}};

      term_2.c--;
      term_2.polynomial = term_2.polynomial * B_00 * term_2.a;
      term_2.a--;

      return {term_1, term_2, term_3};
    }
  }
}

std::vector<IntegralInfo>
horizontal_recursion_relation(const std::vector<IntegralInfo> & info) {

  std::vector<IntegralInfo> result{};

  for (const auto & i_info: info) {
    if (i_info.b > 0) {
      const auto b_horiz = i_info.horizontal_recursion_relation_b_to_a();
      const auto b_recursive = horizontal_recursion_relation(b_horiz);
      result.insert(result.end(), b_recursive.begin(), b_recursive.end());

    } else if (i_info.d > 0) {
      const auto d_horiz = i_info.horizontal_recursion_relation_c_to_d();
      const auto d_recursive = horizontal_recursion_relation(d_horiz);
      result.insert(result.end(), d_recursive.begin(), d_recursive.end());
    } else {
      result.push_back(i_info);
    }
  }

  return result;

}

std::vector<IntegralInfo>
vertical_recursion_relation(const std::vector<IntegralInfo> & info) {

  std::vector<IntegralInfo> result{};

  for (const auto & i_info: info) {
    assert(i_info.b == 0 && i_info.d == 0);
    if (i_info.a > 0) {
      const auto a_vert = i_info.vertical_recursion_relation_a();
      const auto a_recursive = vertical_recursion_relation(a_vert);
      result.insert(result.end(), a_recursive.begin(), a_recursive.end());

    } else if (i_info.c > 0) {
      const auto c_vert = i_info.vertical_recursion_relation_c();
      const auto c_recursive = vertical_recursion_relation(c_vert);
      result.insert(result.end(), c_recursive.begin(), c_recursive.end());
    } else {
      result.push_back(i_info);
    }
  }

  return result;

}


RysPolynomial
reduce_to_rys_polynomial(const IntegralInfo & info) {
  const auto reduced = vertical_recursion_relation(
      horizontal_recursion_relation({info}));

  const int total_angular_momentum = info.a + info.b + info.c + info.d + 1;

  RysPolynomial result = {arma::vec(total_angular_momentum, arma::fill::zeros)};

  for (const auto & i_info: reduced) {
    const arma::uword n_elem = i_info.polynomial.coef.n_elem;
    result.coef(arma::span(0, n_elem - 1)) += i_info.polynomial.coef;
  }

  return result;
}

double electron_repulsive_integral(const ERI & eri_info) {
  const double p = eri_info.A.exponent + eri_info.B.exponent;
  const arma::vec3 P =
      (eri_info.A.exponent * eri_info.A.center +
       eri_info.B.exponent * eri_info.B.center) / p;
  const arma::vec3 from_B_to_A = eri_info.A.center - eri_info.B.center;

  const double exponential_prefactor_AB = std::exp(
      -eri_info.A.exponent * eri_info.B.exponent / p *
      arma::dot(from_B_to_A, from_B_to_A));

  const double q = eri_info.C.exponent + eri_info.D.exponent;
  const arma::vec3 Q =
      (eri_info.C.exponent * eri_info.C.center +
       eri_info.D.exponent * eri_info.D.center) / q;
  const arma::vec3 from_D_to_C = eri_info.C.center - eri_info.D.center;

  const double exponential_prefactor_CD = std::exp(
      -eri_info.C.exponent * eri_info.D.exponent / q *
      arma::dot(from_D_to_C, from_D_to_C));

  const double rho = p * q / (p + q);
  const arma::vec3 from_Q_to_P = P - Q;
  const double T = rho * arma::dot(from_Q_to_P, from_Q_to_P);

  const IntegralInfo I_x{{arma::vec{1.0}},
                         eri_info.A.angular[0],
                         eri_info.B.angular[0],
                         eri_info.C.angular[0],
                         eri_info.D.angular[0],
                         p, P[0], q, Q[0],
                         eri_info.A.center[0],
                         eri_info.B.center[0],
                         eri_info.C.center[0],
                         eri_info.D.center[0],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent};

  const IntegralInfo I_y{{arma::vec{1.0}},
                         eri_info.A.angular[1],
                         eri_info.B.angular[1],
                         eri_info.C.angular[1],
                         eri_info.D.angular[1],
                         p, P[1], q, Q[1],
                         eri_info.A.center[1],
                         eri_info.B.center[1],
                         eri_info.C.center[1],
                         eri_info.D.center[1],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent};


  const double prefactor_for_eri =
      2.0 * std::pow(M_PI, 2.5) / p / q / std::sqrt(p + q) *
      exponential_prefactor_AB * exponential_prefactor_CD;

  const IntegralInfo I_z{{arma::vec{prefactor_for_eri}},
                         eri_info.A.angular[2],
                         eri_info.B.angular[2],
                         eri_info.C.angular[2],
                         eri_info.D.angular[2],
                         p, P[2], q, Q[2],
                         eri_info.A.center[2],
                         eri_info.B.center[2],
                         eri_info.C.center[2],
                         eri_info.D.center[2],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent};

  const RysPolynomial I_z_polynomial = reduce_to_rys_polynomial(I_z);
  const RysPolynomial I_x_polynomial = reduce_to_rys_polynomial(I_x);
  const RysPolynomial I_y_polynomial = reduce_to_rys_polynomial(I_y);

  const RysPolynomial multiplied =
      I_x_polynomial * I_y_polynomial * I_z_polynomial;
  const int n_rys_roots = (multiplied.coef.n_elem + 1) / 2;

  arma::vec u(n_rys_roots);
  arma::vec w(n_rys_roots);

  CINTrys_roots(n_rys_roots, T, u.memptr(), w.memptr());

  const arma::vec t_squared = u / (u + 1.0);

  double sum = 0;
  for (int i = 0; i < n_rys_roots; i++) {
    sum += multiplied(t_squared[i]) * w[i];
  }

  return sum
         * eri_info.A.coef
         * eri_info.B.coef
         * eri_info.C.coef
         * eri_info.D.coef;
}
}

namespace integral::rys_quadrature::gradient {

std::vector<IntegralInfo>
parse_gradient_a(const IntegralInfo & info) {

  auto term_2 = info;
  term_2.grad_a--;
  term_2.a++;
  term_2.polynomial = term_2.polynomial * 2.0 * term_2.alpha;

  if (info.a <= 0) {
    return {term_2};
  } else {
    auto term_1 = info;
    term_1.grad_a--;
    term_1.polynomial = term_1.polynomial * (-info.a);
    term_1.a--;

    return {term_1, term_2};
  }

}

std::vector<IntegralInfo>
parse_gradient_b(const IntegralInfo & info) {

  auto term_2 = info;
  term_2.grad_b--;
  term_2.b++;
  term_2.polynomial = term_2.polynomial * 2.0 * term_2.beta;

  if (info.b <= 0) {
    return {term_2};
  } else {
    auto term_1 = info;
    term_1.grad_b--;
    term_1.polynomial = term_1.polynomial * (-info.b);
    term_1.b--;

    return {term_1, term_2};
  }

}

std::vector<IntegralInfo>
parse_gradient_c(const IntegralInfo & info) {

  auto term_2 = info;
  term_2.grad_c--;
  term_2.c++;
  term_2.polynomial = term_2.polynomial * 2.0 * term_2.gamma;

  if (info.c <= 0) {
    return {term_2};
  } else {
    auto term_1 = info;
    term_1.grad_c--;
    term_1.polynomial = term_1.polynomial * (-info.c);
    term_1.c--;

    return {term_1, term_2};
  }

}

std::vector<IntegralInfo>
parse_gradient_d(const IntegralInfo & info) {

  auto term_2 = info;
  term_2.grad_d--;
  term_2.d++;
  term_2.polynomial = term_2.polynomial * 2.0 * term_2.delta;

  if (info.d <= 0) {
    return {term_2};
  } else {
    auto term_1 = info;
    term_1.grad_d--;
    term_1.polynomial = term_1.polynomial * (-info.d);
    term_1.d--;

    return {term_1, term_2};
  }

}

std::vector<IntegralInfo>
parse_gradient(const std::vector<IntegralInfo> & info) {
  std::vector<IntegralInfo> result{};

  for (const auto & i_info: info) {

    if (i_info.grad_a > 0) {
      const auto after_gradient = parse_gradient_a(i_info);
      const auto after_iteration = parse_gradient(after_gradient);
      result.insert(result.end(), after_iteration.begin(),
                    after_iteration.end());
    } else if (i_info.grad_b > 0) {
      const auto after_gradient = parse_gradient_b(i_info);
      const auto after_iteration = parse_gradient(after_gradient);
      result.insert(result.end(), after_iteration.begin(),
                    after_iteration.end());
    } else if (i_info.grad_c > 0) {
      const auto after_gradient = parse_gradient_c(i_info);
      const auto after_iteration = parse_gradient(after_gradient);
      result.insert(result.end(), after_iteration.begin(),
                    after_iteration.end());
    } else if (i_info.grad_d > 0) {
      const auto after_gradient = parse_gradient_d(i_info);
      const auto after_iteration = parse_gradient(after_gradient);
      result.insert(result.end(), after_iteration.begin(),
                    after_iteration.end());
    } else {
      result.push_back(i_info);
    }
  }

  return result;
}

RysPolynomial
reduce_to_rys_polynomial(const IntegralInfo & info) {

  const auto reduced = vertical_recursion_relation(
      horizontal_recursion_relation(parse_gradient({info})));

  const int total_angular_momentum =
      info.a + info.b + info.c + info.d +
      info.grad_a + info.grad_b + info.grad_c + info.grad_d + 1;

  RysPolynomial result = {arma::vec(total_angular_momentum, arma::fill::zeros)};

  for (const auto & i_info: reduced) {
    const arma::uword n_elem = i_info.polynomial.coef.n_elem;
    result.coef(arma::span(0, n_elem - 1)) += i_info.polynomial.coef;
  }

  return result;
}

double electron_repulsive_integral(const ERI & eri_info,
                                   const arma::Mat<int>::fixed<3, 4> & derivative_operator) {
  const double p = eri_info.A.exponent + eri_info.B.exponent;
  const arma::vec3 P =
      (eri_info.A.exponent * eri_info.A.center +
       eri_info.B.exponent * eri_info.B.center) / p;
  const arma::vec3 from_B_to_A = eri_info.A.center - eri_info.B.center;

  const double exponential_prefactor_AB = std::exp(
      -eri_info.A.exponent * eri_info.B.exponent / p *
      arma::dot(from_B_to_A, from_B_to_A));

  const double q = eri_info.C.exponent + eri_info.D.exponent;
  const arma::vec3 Q =
      (eri_info.C.exponent * eri_info.C.center +
       eri_info.D.exponent * eri_info.D.center) / q;
  const arma::vec3 from_D_to_C = eri_info.C.center - eri_info.D.center;

  const double exponential_prefactor_CD = std::exp(
      -eri_info.C.exponent * eri_info.D.exponent / q *
      arma::dot(from_D_to_C, from_D_to_C));

  const double rho = p * q / (p + q);
  const arma::vec3 from_Q_to_P = P - Q;
  const double T = rho * arma::dot(from_Q_to_P, from_Q_to_P);

  const IntegralInfo I_x{{arma::vec{1.0}},
                         eri_info.A.angular[0],
                         eri_info.B.angular[0],
                         eri_info.C.angular[0],
                         eri_info.D.angular[0],
                         p, P[0], q, Q[0],
                         eri_info.A.center[0],
                         eri_info.B.center[0],
                         eri_info.C.center[0],
                         eri_info.D.center[0],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent,
                         derivative_operator(0, 0),
                         derivative_operator(0, 1),
                         derivative_operator(0, 2),
                         derivative_operator(0, 3)};

  const IntegralInfo I_y{{arma::vec{1.0}},
                         eri_info.A.angular[1],
                         eri_info.B.angular[1],
                         eri_info.C.angular[1],
                         eri_info.D.angular[1],
                         p, P[1], q, Q[1],
                         eri_info.A.center[1],
                         eri_info.B.center[1],
                         eri_info.C.center[1],
                         eri_info.D.center[1],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent,
                         derivative_operator(1, 0),
                         derivative_operator(1, 1),
                         derivative_operator(1, 2),
                         derivative_operator(1, 3)};


  const double prefactor_for_eri =
      2.0 * std::pow(M_PI, 2.5) / p / q / std::sqrt(p + q) *
      exponential_prefactor_AB * exponential_prefactor_CD;

  const IntegralInfo I_z{{arma::vec{prefactor_for_eri}},
                         eri_info.A.angular[2],
                         eri_info.B.angular[2],
                         eri_info.C.angular[2],
                         eri_info.D.angular[2],
                         p, P[2], q, Q[2],
                         eri_info.A.center[2],
                         eri_info.B.center[2],
                         eri_info.C.center[2],
                         eri_info.D.center[2],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent,
                         derivative_operator(2, 0),
                         derivative_operator(2, 1),
                         derivative_operator(2, 2),
                         derivative_operator(2, 3)};

  const RysPolynomial I_z_polynomial = gradient::reduce_to_rys_polynomial(I_z);
  const RysPolynomial I_x_polynomial = gradient::reduce_to_rys_polynomial(I_x);
  const RysPolynomial I_y_polynomial = gradient::reduce_to_rys_polynomial(I_y);

  const RysPolynomial multiplied =
      I_x_polynomial * I_y_polynomial * I_z_polynomial;
  const int n_rys_roots = (multiplied.coef.n_elem + 1) / 2;

  arma::vec u(n_rys_roots);
  arma::vec w(n_rys_roots);

  CINTrys_roots(n_rys_roots, T, u.memptr(), w.memptr());

  const arma::vec t_squared = u / (u + 1.0);

  double sum = 0;
  for (int i = 0; i < n_rys_roots; i++) {
    sum += multiplied(t_squared[i]) * w[i];
  }

  return sum
         * eri_info.A.coef
         * eri_info.B.coef
         * eri_info.C.coef
         * eri_info.D.coef;
}
}

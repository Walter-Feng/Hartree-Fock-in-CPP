#include "integral.h"
#include "../integral.h"

extern "C" {
#include <rys_roots.h>
}

namespace hfincpp::integral::faster_kernel::rys_quadrature {

template<typename T, int a, int b>
inline
void recursive_kernel(double * cache,
                      const T C00,
                      const T D00,
                      const T B01,
                      const T B10,
                      const T B00,
                      const T prefactor = 1) {
  T second_previous, previous = 0, present = prefactor;
  cache[0] = prefactor;
  if constexpr (a > 0) {
#pragma unroll
    for (int a_index = 0; a_index < a; a_index++) {
      second_previous = previous;
      previous = present;
      present = B10 * a_index * second_previous + C00 * previous;
      cache[a_index + 1] = present;
    }
  }

  if constexpr (b > 0) {
    present = prefactor;
#pragma unroll
    for (int b_index = 0; b_index < b; b_index++) {
      second_previous = previous;
      previous = present;
      present = B01 * b_index * second_previous + D00 * previous;
      cache[(a + 1) * (b_index + 1)] = present;

      if constexpr (a > 0) {
        T second_previous_along_a, previous_along_a = 0, present_along_a = present;
#pragma unroll
        for (int a_index = 0; a_index < a; a_index++) {
          second_previous_along_a = previous_along_a;
          previous_along_a = present_along_a;
          present_along_a =
              B10 * a_index * second_previous_along_a + C00 * previous_along_a +
              B00 * (b_index + 1) * cache[a_index + (a + 1) * b_index];
          cache[a_index + 1 + (a + 1) * (b_index + 1)] = present_along_a;
        }
      }
    }
  }
}

constexpr inline size_t binomial_coefficient(size_t n, size_t k) noexcept {
  return
      (k > n) ? 0 :          // out of range
      (k == 0 || k == n) ? 1 :          // edge
      (k == 1 || k == n - 1) ? n :          // first
      binomial_coefficient(n - 1, k - 1) * n / k;   // recursive
}

double binomial_expand(int i,
                       int a,
                       int b,
                       double PA,
                       double PB);

namespace to_conventional {
template<int a, int b>
double electron_repulsive_integral(const ERI & eri_info) {
  const double p = eri_info.A.exponent + eri_info.B.exponent;
  const arma::vec3 P =
      (eri_info.A.exponent * eri_info.A.center +
       eri_info.B.exponent * eri_info.B.center) / p;
  const arma::vec3 from_B_to_A = eri_info.A.center - eri_info.B.center;
  const arma::vec3 from_A_to_P = P - eri_info.A.center;
  const arma::vec3 from_B_to_P = P - eri_info.B.center;

  const double exponential_prefactor_AB = std::exp(
      -eri_info.A.exponent * eri_info.B.exponent / p *
      arma::dot(from_B_to_A, from_B_to_A));

  const double q = eri_info.C.exponent + eri_info.D.exponent;
  const arma::vec3 Q =
      (eri_info.C.exponent * eri_info.C.center +
       eri_info.D.exponent * eri_info.D.center) / q;
  const arma::vec3 from_D_to_C = eri_info.C.center - eri_info.D.center;
  const arma::vec3 from_C_to_Q = Q - eri_info.C.center;
  const arma::vec3 from_D_to_Q = Q - eri_info.D.center;

  const double exponential_prefactor_CD = std::exp(
      -eri_info.C.exponent * eri_info.D.exponent / q *
      arma::dot(from_D_to_C, from_D_to_C));

  const double rho = p * q / (p + q);
  const arma::vec3 from_Q_to_P = P - Q;
  const double T = rho * arma::dot(from_Q_to_P, from_Q_to_P);

  const double prefactor_for_eri =
      2.0 * std::pow(M_PI, 2.5) / p / q / std::sqrt(p + q) *
      exponential_prefactor_AB * exponential_prefactor_CD
      * eri_info.A.coef
      * eri_info.B.coef
      * eri_info.C.coef
      * eri_info.D.coef;

  assert(
      arma::sum(eri_info.A.angular) + arma::sum(eri_info.B.angular) == a);
  assert(
      arma::sum(eri_info.C.angular) + arma::sum(eri_info.D.angular) == b);

  const int n_rys_roots = (a + b) / 2 + 1;

  arma::vec u(n_rys_roots);
  arma::vec w(n_rys_roots);

  CINTrys_roots(n_rys_roots, T, u.memptr(), w.memptr());

  const arma::vec t_squared_list = u / (u + 1.0);

  const arma::vec exponent_weighted_t_squared = t_squared_list / (p + q);

  double C00, D00, B01, B10, B00, PA, PB, QC, QD;
  arma::mat integral_table(a + 1, b + 1);

  double rys_sum = 0;
  for (int i_rys = 0; i_rys < n_rys_roots; i_rys++) {
    double t_squared = exponent_weighted_t_squared(i_rys);
    double weight = w(i_rys);

    double integral_per_grid = weight;

    B00 = t_squared / 2.0;
    B10 = 0.5 / p - B00 * q / p;
    B01 = 0.5 / q - B00 * p / q;

    for (int xyz_index = 0; xyz_index < 3; xyz_index++) {
      C00 = -t_squared * q * from_Q_to_P(xyz_index);
      D00 = t_squared * p * from_Q_to_P(xyz_index);

      PA = from_A_to_P(xyz_index);
      PB = from_B_to_P(xyz_index);
      QC = from_C_to_Q(xyz_index);
      QD = from_D_to_Q(xyz_index);

      double xyz_component = 0;
      recursive_kernel<double, a, b>(integral_table.memptr(),
                                     C00, D00, B01, B10, B00);

      int i = eri_info.A.angular(xyz_index);
      int j = eri_info.B.angular(xyz_index);
      int k = eri_info.C.angular(xyz_index);
      int l = eri_info.D.angular(xyz_index);

      for (int ij = 0; ij <= i + j; ij++) {
        for (int kl = 0; kl <= k + l; kl++) {
          xyz_component += integral_table(ij, kl)
                           * binomial_expand(ij, i, j, PA, PB)
                           * binomial_expand(kl, k, l, QC, QD);

        }
      }

      integral_per_grid *= xyz_component;
    }

    rys_sum += integral_per_grid;
  }
  return rys_sum * prefactor_for_eri;
}
}

}


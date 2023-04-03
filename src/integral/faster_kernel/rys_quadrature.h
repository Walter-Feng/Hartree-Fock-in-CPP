#include "integral.h"
#include "../integral.h"

extern "C" {
#include <rys_roots.h>
}

namespace hfincpp::integral::faster_kernel::rys_quadrature {

template<int a, int b, int stride>
void electron_repulsive_integral_kernel(double * cache,
                                        const double C00,
                                        const double D00,
                                        const double B01,
                                        const double B10,
                                        const double B00) {
  int one = 1;
  int two = 2;
  int a_width = a + one;
  int minus_stride = -stride;
  int minus_b_stride = a_width * minus_stride;
  *cache = 1.;
  cache += stride;

  for (int i = one; i <= a; i++) {
    *cache =
        C00 * cache[minus_stride] + B10 * (i - one) * cache[two * minus_stride];

    cache += stride;
  }
  for (int j = one; j <= b; j++) {
    for (int i = 0; i <= a; i++) {
      *cache =
          D00 * cache[minus_b_stride]
          + B01 * (j-one) * cache[two * minus_b_stride]
          + B00 * i * cache[minus_stride + minus_b_stride];

      cache += stride;
    }
  }
}

int binomial_coefficient(int n, int r);

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

  assert(arma::sum(eri_info.A.angular) + arma::sum(eri_info.B.angular) == a);
  assert(arma::sum(eri_info.C.angular) + arma::sum(eri_info.D.angular) == b);

  const int n_rys_roots = (a + b + 1) / 2;

  arma::vec u(n_rys_roots);
  arma::vec w(n_rys_roots);

  CINTrys_roots(n_rys_roots, T, u.memptr(), w.memptr());

  const arma::vec t_squared_list = u / (u + 1.0);

  const arma::vec exponent_weighted_t_squared = t_squared_list / (p+q);

  double C00, D00, B01, B10, B00, PA, PB, QC, QD;
  arma::mat integral_table(a + 1, b + 1);

  double rys_sum = 0;
  for(int i_rys=0; i_rys<n_rys_roots; i_rys++) {
    double t_squared = exponent_weighted_t_squared(i_rys);
    double weight = w(i_rys);

    double integral_per_grid = weight;

    B00 = t_squared / 2.0;
    B10 = 0.5 / p - B00 * q / p;
    B01 = 0.5 / q - B00 * p / q;

    for(int xyz_index=0; xyz_index<3; xyz_index++) {
      C00 = -  t_squared * q * from_Q_to_P(xyz_index);
      D00 = t_squared * p * from_Q_to_P(xyz_index);

      PA = from_A_to_P(xyz_index);
      PB = from_B_to_P(xyz_index);
      QC = from_C_to_Q(xyz_index);
      QD = from_D_to_Q(xyz_index);

      double xyz_component = 0;
      electron_repulsive_integral_kernel<a,b,1>(integral_table.memptr(),
                                                C00, D00, B01, B10, B00);

      int i = eri_info.A.angular(xyz_index);
      int j = eri_info.B.angular(xyz_index);
      int k = eri_info.C.angular(xyz_index);
      int l = eri_info.D.angular(xyz_index);

      for(int ij=0; ij<=i+j; ij++) {
        for(int kl=0; kl<=k+l; kl++) {
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


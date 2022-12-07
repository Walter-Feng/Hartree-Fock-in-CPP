#include "obara_saika.h"

#include <boost/math/special_functions.hpp>
#include <cmath>

namespace hfincpp::integral::obara_saika {

double Boys(double x, int m) {
  if (std::abs(x) < 1e-8) return (1.0 / (1.0 + 2.0 * m));
  else
    return 0.5 * std::pow(x, -0.5 - m) *
           (boost::math::tgamma(0.5 + m) - boost::math::tgamma(0.5 + m, x));
}

/*
 * The functions below are directly based on a thesis
 * May, Andrew James. Density fitting in explicitly correlated
 * electronic structure theory. Diss. University of Bristol, 2006.
 * They can also be referenced from the paper
 * Obara, Shigeru, and A. Saika. "Efficient recursive computation of
 * molecular integrals over Cartesian Gaussian functions."
 * The Journal of chemical physics 84.7 (1986): 3963-3974.
 * Only the overlap functions are used in actual code.
 */

//Calculate the Overlap Integral of two gaussian function
double overlap_integral(const double ra[3],
                        const double rb[3],
                        const int ax, const int ay, const int az,
                        const int bx, const int by, const int bz,
                        const double alpha, const double beta) {

  // if one of the angular number goes below zero,
  // it means it will not have contribution - the same as giving derivation to a constant
  if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0) return 0;

    //Provide recurrence relation
  else if (ax > 0)
    return ((alpha * ra[0] + beta * rb[0]) / (alpha + beta) - ra[0]) *
           overlap_integral(ra, rb, ax - 1, ay, az, bx, by, bz, alpha, beta) +
           (ax - 1) / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax - 2, ay, az, bx, by, bz, alpha, beta) +
           bx / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha,
                            beta);

  else if (ay > 0)
    return ((alpha * ra[1] + beta * rb[1]) / (alpha + beta) - ra[1]) *
           overlap_integral(ra, rb, ax, ay - 1, az, bx, by, bz, alpha, beta) +
           (ay - 1) / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay - 2, az, bx, by, bz, alpha, beta) +
           by / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha,
                            beta);

  else if (az > 0)
    return ((alpha * ra[2] + beta * rb[2]) / (alpha + beta) - ra[2]) *
           overlap_integral(ra, rb, ax, ay, az - 1, bx, by, bz, alpha, beta) +
           (az - 1) / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay, az - 2, bx, by, bz, alpha, beta) +
           bz / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha,
                            beta);

  else if (bx > 0)
    return ((alpha * ra[0] + beta * rb[0]) / (alpha + beta) - rb[0]) *
           overlap_integral(ra, rb, ax, ay, az, bx - 1, by, bz, alpha, beta) +
           ax / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha,
                            beta) +
           (bx - 1) / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay, az, bx - 2, by, bz, alpha, beta);

  else if (by > 0)
    return ((alpha * ra[1] + beta * rb[1]) / (alpha + beta) - rb[1]) *
           overlap_integral(ra, rb, ax, ay, az, bx, by - 1, bz, alpha, beta) +
           ay / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha,
                            beta) +
           (by - 1) / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay, az, bx, by - 2, bz, alpha, beta);

  else if (bz > 0)
    return ((alpha * ra[2] + beta * rb[2]) / (alpha + beta) - rb[2]) *
           overlap_integral(ra, rb, ax, ay, az, bx, by, bz - 1, alpha, beta) +
           az / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha,
                            beta) +
           (bz - 1) / 2.0 / (alpha + beta) *
           overlap_integral(ra, rb, ax, ay, az, bx, by, bz - 2, alpha, beta);

    //giving the starting point
  else
    return sqrt(M_PI / (alpha + beta)) * M_PI / (alpha + beta) *
           exp(-alpha * beta / (alpha + beta) *
               (pow(ra[0] - rb[0], 2) + pow(ra[1] - rb[1], 2) +
                pow(ra[2] - rb[2], 2)));
}

double
overlap_integral(const GaussianFunction & A, const GaussianFunction & B) {
  return overlap_integral(A.center.memptr(), B.center.memptr(),
                          A.angular[0], A.angular[1], A.angular[2],
                          B.angular[0], B.angular[1], B.angular[2],
                          A.exponent, B.exponent) * A.coef * B.coef;
}

arma::mat overlap_integral(const basis::Basis & basis) {

  arma::mat overlap(basis.n_functions(), basis.n_functions());

#pragma omp parallel for collapse(2)
  for (int i = 0; i < basis.n_functions(); i++) {
    for (int j = i; j < basis.n_functions(); j++) {
      const auto & function_i = basis.functions[i];
      const auto n_gto_from_i = function_i.coefficients.n_elem;
      const auto & function_j = basis.functions[j];
      const auto n_gto_from_j = function_j.coefficients.n_elem;

      double value = 0;
      for (arma::uword gto_i = 0; gto_i < n_gto_from_i; gto_i++) {
        for (arma::uword gto_j = 0; gto_j < n_gto_from_j; gto_j++) {
          const GaussianFunction gto_function_i{function_i.center,
                                                function_i.angular,
                                                function_i.exponents(gto_i),
                                                function_i.coefficients(gto_i)};

          const GaussianFunction gto_function_j{function_j.center,
                                                function_j.angular,
                                                function_j.exponents(gto_j),
                                                function_j.coefficients(gto_j)};

          value += overlap_integral(gto_function_i, gto_function_j);
        }
      }

      overlap(i, j) = value;
      overlap(j, i) = value;
    }
  }

  return overlap;
}

arma::mat kinetic_integral(const basis::Basis & basis) {

  arma::mat kinetic(basis.n_functions(), basis.n_functions());

#pragma omp parallel for collapse(2)
  for (int i = 0; i < basis.n_functions(); i++) {
    for (int j = i; j < basis.n_functions(); j++) {
      const auto & function_i = basis.functions[i];
      const auto n_gto_from_i = function_i.coefficients.n_elem;
      const auto & function_j = basis.functions[j];
      const auto n_gto_from_j = function_j.coefficients.n_elem;

      double value = 0;
      for (arma::uword gto_i = 0; gto_i < n_gto_from_i; gto_i++) {
        for (arma::uword gto_j = 0; gto_j < n_gto_from_j; gto_j++) {
          const GaussianFunction gto_function_i{function_i.center,
                                                function_i.angular,
                                                function_i.exponents(gto_i),
                                                function_i.coefficients(gto_i)};

          const GaussianFunction gto_function_j{function_j.center,
                                                function_j.angular,
                                                function_j.exponents(gto_j),
                                                function_j.coefficients(gto_j)};

          const auto laplace_operator_on_j = gto_function_j.laplace();

          for (const auto & gto_k: laplace_operator_on_j) {
            value += overlap_integral(gto_function_i, gto_k);
          }
        }
      }

      kinetic(i, j) = -0.5 * value;
      kinetic(j, i) = -0.5 * value;
    }
  }

  return kinetic;
}

struct ERIIntermediate {
  double zeta;
  double eta;
  double chi;
  double rho;
  arma::vec3 P;
  arma::vec3 Q;
  arma::vec3 W;
};

double
electron_repulsive_integral(const GaussianFunction & A,
                            const GaussianFunction & B,
                            const GaussianFunction & C,
                            const GaussianFunction & D,
                            const ERIIntermediate & intermediate,
                            const int m) {
  if (arma::any(A.angular < 0)
      || arma::any(B.angular < 0)
      || arma::any(C.angular < 0)
      || arma::any(D.angular < 0)) {
    return 0;
  }

  const auto & P = intermediate.P;
  const auto & Q = intermediate.Q;
  const auto & W = intermediate.W;
  const auto & zeta = intermediate.zeta;
  const auto & eta = intermediate.eta;
  const auto & chi = intermediate.chi;
  const auto & rho = intermediate.rho;


  for (int i = 0; i < 3; i++) {
    auto reduced_A = A;
    auto reduced_B = B;
    auto reduced_C = C;
    auto reduced_D = D;

    reduced_A.angular[i]--;
    reduced_B.angular[i]--;
    reduced_C.angular[i]--;
    reduced_D.angular[i]--;

    if (A.angular[i] > 0) {
      auto double_reduced_A = reduced_A;
      double_reduced_A.angular[i]--;

      return (P[i] - A.center[i]) *
             electron_repulsive_integral(reduced_A, B, C, D, intermediate, m)
             + (W[i] - P[i]) *
               electron_repulsive_integral(reduced_A, B, C, D, intermediate,
                                           m + 1)

             + 0.5 / zeta * reduced_A.angular[i] * (
          electron_repulsive_integral(double_reduced_A, B, C, D, intermediate,
                                      m)
          - rho / zeta *
            electron_repulsive_integral(double_reduced_A, B, C, D, intermediate,
                                        m + 1)
      ) + 0.5 / zeta * B.angular[i] * (
          electron_repulsive_integral(reduced_A, reduced_B, C, D, intermediate,
                                      m)
          - rho / zeta *
            electron_repulsive_integral(reduced_A, reduced_B, C, D,
                                        intermediate, m + 1)
      ) + 0.5 / chi * C.angular[i] *
          electron_repulsive_integral(reduced_A, B, reduced_C, D, intermediate,
                                      m + 1)

             + 0.5 / chi * D.angular[i] *
               electron_repulsive_integral(reduced_A, B, C, reduced_D,
                                           intermediate,
                                           m + 1);
    }

    if (B.angular[i] > 0) {
      auto double_reduced_B = reduced_B;
      double_reduced_B.angular[i]--;

      return (P[i] - B.center[i]) *
             electron_repulsive_integral(A, reduced_B, C, D, intermediate, m)
             + (W[i] - P[i]) *
               electron_repulsive_integral(A, reduced_B, C, D, intermediate,
                                           m + 1)

             + 0.5 / zeta * reduced_B.angular[i] * (
          electron_repulsive_integral(A, double_reduced_B, C, D, intermediate,
                                      m)
          - rho / zeta *
            electron_repulsive_integral(A, double_reduced_B, C, D, intermediate,
                                        m + 1)
      ) + 0.5 / zeta * A.angular[i] * (
          electron_repulsive_integral(reduced_A, reduced_B, C, D, intermediate,
                                      m)
          - rho / zeta *
            electron_repulsive_integral(reduced_A, reduced_B, C, D,
                                        intermediate, m + 1)
      ) + 0.5 / chi * C.angular[i] *
          electron_repulsive_integral(A, reduced_B, reduced_C, D, intermediate,
                                      m + 1)

             + 0.5 / chi * D.angular[i] *
               electron_repulsive_integral(A, reduced_B, C, reduced_D,
                                           intermediate,
                                           m + 1);
    }

    if (C.angular[i] > 0) {
      auto double_reduced_C = reduced_C;
      double_reduced_C.angular[i]--;

      return (Q[i] - C.center[i]) *
             electron_repulsive_integral(A, B, reduced_C, D, intermediate, m)
             + (W[i] - Q[i]) *
               electron_repulsive_integral(A, B, reduced_C, D, intermediate,
                                           m + 1)

             + 0.5 / eta * reduced_C.angular[i] * (
          electron_repulsive_integral(A, B, double_reduced_C, D, intermediate,
                                      m)
          - rho / eta *
            electron_repulsive_integral(A, B, double_reduced_C, D, intermediate,
                                        m + 1)
      ) + 0.5 / eta * D.angular[i] * (
          electron_repulsive_integral(A, B, reduced_C, reduced_D, intermediate,
                                      m)
          - rho / eta *
            electron_repulsive_integral(A, B, reduced_C, reduced_D,
                                        intermediate, m + 1)
      ) + 0.5 / chi * A.angular[i] *
          electron_repulsive_integral(reduced_A, B, reduced_C, D, intermediate,
                                      m + 1)

             + 0.5 / chi * B.angular[i] *
               electron_repulsive_integral(A, reduced_B, reduced_C, D,
                                           intermediate,
                                           m + 1);
    }

    if (D.angular[i] > 0) {
      auto double_reduced_D = reduced_D;
      double_reduced_D.angular[i]--;

      return (Q[i] - D.center[i]) *
             electron_repulsive_integral(A, B, C, reduced_D, intermediate, m)
             + (W[i] - Q[i]) *
               electron_repulsive_integral(A, B, C, reduced_D, intermediate,
                                           m + 1)

             + 0.5 / eta * reduced_D.angular[i] * (
          electron_repulsive_integral(A, B, C, double_reduced_D, intermediate,
                                      m)
          - rho / eta *
            electron_repulsive_integral(A, B, C, double_reduced_D, intermediate,
                                        m + 1)
      ) + 0.5 / eta * C.angular[i] * (
          electron_repulsive_integral(A, B, reduced_C, reduced_D, intermediate,
                                      m)
          - rho / eta *
            electron_repulsive_integral(A, B, reduced_C, reduced_D,
                                        intermediate, m + 1)
      ) + 0.5 / chi * A.angular[i] *
          electron_repulsive_integral(reduced_A, B, C, reduced_D, intermediate,
                                      m + 1)

             + 0.5 / chi * B.angular[i] *
               electron_repulsive_integral(A, reduced_B, C, reduced_D,
                                           intermediate,
                                           m + 1);
    }
  }

  const arma::vec3 from_B_to_A = A.center - B.center;
  const arma::vec3 from_D_to_C = C.center - D.center;
  const arma::vec3 from_Q_to_P = P - Q;
  const double AB_norm_squared = arma::sum(arma::square(from_B_to_A));
  const double CD_norm_squared = arma::sum(arma::square(from_D_to_C));
  const double PQ_norm_squared = arma::sum(arma::square(from_Q_to_P));
  const double xi = A.exponent * B.exponent / (A.exponent + B.exponent);
  const double xi_prime = C.exponent * D.exponent / (C.exponent + D.exponent);

  assert((arma::any(A.angular == 0)
          && arma::any(B.angular == 0)
          && arma::any(C.angular == 0)
          && arma::any(D.angular == 0)));
  return 2.0 * std::pow(M_PI, 2.5) / (zeta * eta) / std::sqrt(zeta + eta)
         * std::exp(-xi * AB_norm_squared) *
         std::exp(-xi_prime * CD_norm_squared)
         * Boys(rho * PQ_norm_squared, m);
}

double
electron_repulsive_integral(const GaussianFunction & A,
                            const GaussianFunction & B,
                            const GaussianFunction & C,
                            const GaussianFunction & D) {
  const double zeta = A.exponent + B.exponent;
  const double eta = C.exponent + D.exponent;
  const double chi = eta + zeta;
  const double rho = zeta * eta / chi;
  const arma::vec3 P = (A.exponent * A.center + B.exponent * B.center) /
                       (A.exponent + B.exponent);
  const arma::vec3 Q = (C.exponent * C.center + D.exponent * D.center) /
                       (C.exponent + D.exponent);
  const arma::vec3 W = (zeta * P + eta * Q) / (zeta + eta);

  const auto intermediate = ERIIntermediate{zeta, eta, chi, rho, P, Q, W};
  return electron_repulsive_integral(A, B, C, D, intermediate, 0)
         * A.coef * B.coef * C.coef * D.coef;
}

double
electron_repulsive_integral(const ERI & eri_info) {
  return electron_repulsive_integral(eri_info.A, eri_info.B, eri_info.C,
                                     eri_info.D);
}

arma::mat electron_repulsive_integral(const basis::Basis & basis) {

  const auto n_ao = basis.n_functions();
  const auto pair_size = n_ao * n_ao;

  arma::mat eri(pair_size, pair_size);

  /* omp command here is just an automatic parallelization parsing of for loops
   * which has no technical significance. Forget it.
   * Instead of iterating over n_ao for all (ijkl) indices,
   * we can utilize the symmetry of [ij | kl],
   * [ij|kl] = [ji|kl] = [ij|lk] = [kl|ij]
   * These three "=" represent 2^3=8-fold symmetry.
   * the symmetry [ij|kl] =[kl|ij] is slightly different from having symmetric matrix,
   * as we will also be having symmetry for i <-> j and k <-> l which
   * destroys it.
   * To implement this symmetry we usually use shell-pairs that first generate
   * lower-triangular (ij) as pair, and from this list of pairs we pick up
   * (ij) - pair and (kl) - pair as arguments for the ERI, but in a
   * lower-triangular manner. This will require some indexing that, well,
   * not hard but slightly complicated to implement here.
   * So only i<->j symmetry and k<->l symmetry implemented, 4-fold in total.
   */
#pragma omp parallel for collapse(4)
  for (int i = 0; i < n_ao; i++) {
    for (int j = 0; j <= i; j++) {
      for (int k = 0; k < n_ao; k++) {
        for (int l = 0; l <= k; l++) {
          const auto & function_i = basis.functions[i];
          const auto n_gto_from_i = function_i.coefficients.n_elem;
          const auto & function_j = basis.functions[j];
          const auto n_gto_from_j = function_j.coefficients.n_elem;
          const auto & function_k = basis.functions[k];
          const auto n_gto_from_k = function_k.coefficients.n_elem;
          const auto & function_l = basis.functions[l];
          const auto n_gto_from_l = function_l.coefficients.n_elem;


          double value = 0;
          for (arma::uword gto_i = 0; gto_i < n_gto_from_i; gto_i++) {
            for (arma::uword gto_j = 0; gto_j < n_gto_from_j; gto_j++) {
              for (arma::uword gto_k = 0; gto_k < n_gto_from_k; gto_k++) {
                for (arma::uword gto_l = 0; gto_l < n_gto_from_l; gto_l++) {
                  const GaussianFunction gto_function_i{function_i.center,
                                                        function_i.angular,
                                                        function_i.exponents(
                                                            gto_i),
                                                        function_i.coefficients(
                                                            gto_i)};

                  const GaussianFunction gto_function_j{function_j.center,
                                                        function_j.angular,
                                                        function_j.exponents(
                                                            gto_j),
                                                        function_j.coefficients(
                                                            gto_j)};

                  const GaussianFunction gto_function_k{function_k.center,
                                                        function_k.angular,
                                                        function_k.exponents(
                                                            gto_k),
                                                        function_k.coefficients(
                                                            gto_k)};

                  const GaussianFunction gto_function_l{function_l.center,
                                                        function_l.angular,
                                                        function_l.exponents(
                                                            gto_l),
                                                        function_l.coefficients(
                                                            gto_l)};

                  const ERI eri_info{gto_function_i, gto_function_j,
                                     gto_function_k, gto_function_l};

                  value += electron_repulsive_integral(eri_info);
                }
              }
            }
          }

          eri(i + n_ao * j, k + n_ao * l) = value;
          eri(j + n_ao * i, k + n_ao * l) = value;
          eri(i + n_ao * j, l + n_ao * k) = value;
          eri(j + n_ao * i, l + n_ao * k) = value;
        }
      }
    }
  }

  return eri;

}

// The nuclear attraction integral
double
nuclear_attraction_integral(double ra[3], double rb[3], double rz[3], int ax,
                            int ay, int az,
                            int bx, int by, int bz, double alpha, double beta,
                            int m) {
  double zeta = alpha + beta;
  double P[3];
  double AB, PC;

  int i;

  for (i = 0; i < 3; i++)
    P[i] = (alpha * ra[i] + beta * rb[i]) / zeta;

  AB = 0;
  for (i = 0; i < 3; i++)
    AB += (ra[i] - rb[i]) * (ra[i] - rb[i]);

  PC = 0;
  for (i = 0; i < 3; i++)
    PC += (P[i] - rz[i]) * (P[i] - rz[i]);

  if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0) return 0;

  else if (ax > 0)
    return beta / (alpha + beta) * (rb[0] - ra[0]) *
           nuclear_attraction_integral(ra, rb, rz, ax - 1, ay, az, bx, by, bz,
                                       alpha, beta, m) +
           (rz[0] - ra[0] - beta / (alpha + beta) * (rb[0] - ra[0])) *
           nuclear_attraction_integral(ra, rb, rz, ax - 1, ay, az, bx, by, bz,
                                       alpha, beta,
                                       m + 1) +
           (double) (ax - 1) / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax - 2, ay, az, bx, by, bz,
                                        alpha, beta, m) -
            nuclear_attraction_integral(ra, rb, rz, ax - 2, ay, az, bx, by, bz,
                                        alpha, beta, m + 1)) +
           (double) bx / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax - 1, ay, az, bx - 1, by,
                                        bz, alpha, beta,
                                        m) -
            nuclear_attraction_integral(ra, rb, rz, ax - 1, ay, az, bx - 1, by,
                                        bz,
                                        alpha, beta, m + 1));

  else if (ay > 0)
    return beta / (alpha + beta) * (rb[1] - ra[1]) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay - 1, az, bx, by, bz,
                                       alpha, beta, m) +
           (rz[1] - ra[1] - beta / (alpha + beta) * (rb[1] - ra[1])) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay - 1, az, bx, by, bz,
                                       alpha, beta,
                                       m + 1) +
           (double) (ay - 1) / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay - 2, az, bx, by, bz,
                                        alpha, beta, m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay - 2, az, bx, by, bz,
                                        alpha, beta, m + 1)) +
           (double) by / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay - 1, az, bx, by - 1,
                                        bz, alpha, beta,
                                        m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay - 1, az, bx, by - 1,
                                        bz,
                                        alpha, beta, m + 1));

  else if (az > 0)
    return beta / (alpha + beta) * (rb[2] - ra[2]) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 1, bx, by, bz,
                                       alpha, beta, m) +
           (rz[2] - ra[2] - beta / (alpha + beta) * (rb[2] - ra[2])) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 1, bx, by, bz,
                                       alpha, beta,
                                       m + 1) +
           (double) (az - 1) / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 2, bx, by, bz,
                                        alpha, beta, m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 2, bx, by, bz,
                                        alpha, beta, m + 1)) +
           (double) bz / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 1, bx, by,
                                        bz - 1, alpha, beta,
                                        m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 1, bx, by,
                                        bz - 1,
                                        alpha, beta, m + 1));

  else if (bx > 0)
    return alpha / (alpha + beta) * (ra[0] - rb[0]) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx - 1, by, bz,
                                       alpha, beta, m) +
           (rz[0] - rb[0] - alpha / (alpha + beta) * (ra[0] - rb[0])) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx - 1, by, bz,
                                       alpha, beta,
                                       m + 1) +
           (double) (bx - 1) / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx - 2, by, bz,
                                        alpha, beta, m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx - 2, by, bz,
                                        alpha, beta, m + 1)) +
           (double) ax / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax - 1, ay, az, bx - 1, by,
                                        bz, alpha, beta,
                                        m) -
            nuclear_attraction_integral(ra, rb, rz, ax - 1, ay, az, bx - 1, by,
                                        bz,
                                        alpha, beta, m + 1));

  else if (by > 0)
    return alpha / (alpha + beta) * (ra[1] - rb[1]) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by - 1, bz,
                                       alpha, beta, m) +
           (rz[1] - rb[1] - alpha / (alpha + beta) * (ra[1] - rb[1])) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by - 1, bz,
                                       alpha, beta,
                                       m + 1) +
           (double) (by - 1) / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by - 2, bz,
                                        alpha, beta, m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by - 2, bz,
                                        alpha, beta, m + 1)) +
           (double) ay / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay - 1, az, bx, by - 1,
                                        bz, alpha, beta,
                                        m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay - 1, az, bx, by - 1,
                                        bz,
                                        alpha, beta, m + 1));

  else if (bz > 0)
    return alpha / (alpha + beta) * (ra[2] - rb[2]) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by, bz - 1,
                                       alpha, beta, m) +
           (rz[2] - rb[2] - alpha / (alpha + beta) * (ra[2] - rb[2])) *
           nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by, bz - 1,
                                       alpha, beta,
                                       m + 1) +
           (double) (bz - 1) / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by, bz - 2,
                                        alpha, beta, m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay, az, bx, by, bz - 2,
                                        alpha, beta, m + 1)) +
           (double) az / 2.0 / (alpha + beta) *
           (nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 1, bx, by,
                                        bz - 1, alpha, beta,
                                        m) -
            nuclear_attraction_integral(ra, rb, rz, ax, ay, az - 1, bx, by,
                                        bz - 1,
                                        alpha, beta, m + 1));

  else
    return 2 * M_PI / (alpha + beta) * exp(-alpha * beta / zeta * AB) *
           Boys((alpha + beta) * PC, m);
}

namespace gradient {
arma::cube overlap_integral(const basis::Basis & basis) {

  const arma::uword n_atoms = basis.n_atoms();

  arma::cube overlap(basis.n_functions(),
                     basis.n_functions(),
                     n_atoms * 3,
                     arma::fill::zeros);

  const auto on_atoms = basis.on_atoms();
  for (arma::uword i_atom = 0; i_atom < n_atoms; i_atom++) {
    const auto & on_atom = on_atoms[i_atom];
//#pragma omp parallel for collapse(2)
    for (arma::uword i_function_index = 0;
         i_function_index < on_atom.n_elem; i_function_index++) {
      for (int j = 0; j < basis.n_functions(); j++) {
        if(basis.atom_indices(j) == i_atom) {
          continue;
        }
        const arma::uword i = on_atom(i_function_index);
        const auto & function_i = basis.functions[i];
        const auto n_gto_from_i = function_i.coefficients.n_elem;
        const auto & function_j = basis.functions[j];
        const auto n_gto_from_j = function_j.coefficients.n_elem;

        for (arma::uword gto_i = 0; gto_i < n_gto_from_i; gto_i++) {
          for (arma::uword gto_j = 0; gto_j < n_gto_from_j; gto_j++) {
            for (int xyz_index = 0; xyz_index < 3; xyz_index++) {
              double value = 0;
              const GaussianFunction gto_function_i{function_i.center,
                                                    function_i.angular,
                                                    function_i.exponents(gto_i),
                                                    function_i.coefficients(
                                                        gto_i)};

              const GaussianFunction gto_function_j{function_j.center,
                                                    function_j.angular,
                                                    function_j.exponents(gto_j),
                                                    function_j.coefficients(
                                                        gto_j)};

              const auto gradient_on_i = gto_function_i.gradient(xyz_index);

              for (const auto & gto_k: gradient_on_i) {
                value += obara_saika::overlap_integral(gto_k, gto_function_j);
              }

              overlap(i, j, i_atom * 3 + xyz_index) -= value;
              overlap(j, i, i_atom * 3 + xyz_index) -= value;
            }
          }
        }
      }
    }
  }

  return overlap;
}
}
}
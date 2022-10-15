#include "obara_saika.h"

#include <gsl/gsl_sf.h>
#include <cmath>

namespace integral::obara_saika {

// Calculate the complete Gamma function
double Gamma(double z) {

  const int a = 12;
  static double c_space[12];
  static double * c = nullptr;
  int k;
  double accm;

  if (c == nullptr) {
    double k1_factrl = 1.0; /* (k - 1)!*(-1)^k with 0!==1*/
    c = c_space;
    c[0] = sqrt(2.0 * M_PI);
    for (k = 1; k < a; k++) {
      c[k] = exp(a - k) * pow(a - k, k - 0.5) / k1_factrl;
      k1_factrl *= -k;
    }
  }
  accm = c[0];
  for (k = 1; k < a; k++) {
    accm += c[k] / (z + k);
  }
  accm *= exp(-(z + a)) * pow(z + a, z + 0.5); /* Gamma(z+1) */
  return accm / z;
}

double Boys(double x, int n) {
  if (std::abs(x) < 1e-8) return (1.0 / (1.0 + 2.0 * (double) n));
  else
    return 0.5 * pow(x, -0.5 - n) *
           (Gamma(0.5 + n) - gsl_sf_gamma_inc(0.5 + n, x));
}

double Binomials(int n, int k) {
  int i;
  double result = 1.0;
  if (n == k) return 1;
  if (k == 0) return 1;
  if (k < n) return 0;
  if (k < 0) return 0;
  for (i = 1; i < k; i++)
    result *= (double) n - i;
  for (i = 1; i <= k; i++)
    result = result / (double) i;

  return result;
}

// A special function set for calculating the transformation coefficient
double f(int k, int a, int b, double PA, double PB) {
  int i;
  double result;
  result = 0;

  for (i = 0; i <= k; i++) {
    if (a < i || b - k + i < 0) {
      continue;
    }
    if (abs(PA) < 1e-16 && a - i == 0 && abs(PB) < 1e-19 && b - k + i == 0) {
      result += Binomials(a, i) * Binomials(b, k - i);
    } else if (abs(PA) < 1e-16 && a - i == 0) {
      result += Binomials(a, i) * Binomials(b, k - i) * pow(PB, b - k + i);
    } else if (abs(PB) < 1e-19 && b - k + i == 0) {
      result += Binomials(a, i) * Binomials(a, k - i) * pow(PA, a - i);
    } else {
      result += Binomials(a, i) * Binomials(a, k - i)
                * pow(PA, a - i) * pow(PB, b - k + i);
    }
  }

  return result;
}

std::vector<GaussianFunction>
expand_function_pair(const GaussianFunctionPair & pair) {

  const auto & A = pair.first;
  const auto & B = pair.second;

  const double zeta = A.exponent + B.exponent;
  const double xi = A.exponent * B.exponent / zeta;
  const arma::vec3 P = (A.exponent * A.center + B.exponent * B.center) / zeta;
  const arma::vec3 from_B_to_A = A.center - B.center;
  const double AB_norm_squared = arma::dot(from_B_to_A, from_B_to_A);

  const double exponential_factor = std::exp(-xi * AB_norm_squared);

  const arma::vec3 from_A_to_P = P - A.center;
  const arma::vec3 from_B_to_P = P - B.center;

  std::vector<GaussianFunction> result;
  for (int i = 0; i <= A.angular[0] + B.angular[0]; i++) {
    const double f_x = f(i, A.angular[0], B.angular[0], from_A_to_P[0],
                         from_B_to_P[0]);
    for (int j = 0; j <= A.angular[1] + B.angular[1]; j++) {
      const double f_y = f(j, A.angular[1], B.angular[1], from_A_to_P[1],
                           from_B_to_P[1]);
      for (int k = 0; k <= A.angular[2] + B.angular[2]; k++) {
        const double f_z = f(k, A.angular[2], B.angular[2], from_A_to_P[2],
                             from_B_to_P[2]);

        const double new_coef =
            A.coef * B.coef * exponential_factor * f_x * f_y * f_z;

        result.push_back({P, {i, j, k}, zeta, new_coef});
      }
    }
  }

  return result;
}

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
                          A.exponent, B.exponent);
}

double overlap_integral(const GaussianFunctionPair & pair) {
  return overlap_integral(pair.first, pair.second);
}

//Calculate the Coulomb Integral of two gaussian function
double
electron_repulsive_integral(const double ra[3],
                            const double rb[3],
                            const int ax, const int ay, const int az,
                            const int bx, const int by, const int bz,
                            const double alpha,
                            const double beta,
                            const int m) {
  double zeta = alpha + beta;
  double xi = alpha * beta / zeta;
  double P[3];
  double AB;

  int i;

  for (i = 0; i < 3; i++)
    P[i] = (alpha * ra[i] + beta * rb[i]) / zeta;

  AB = 0;
  for (i = 0; i < 3; i++)
    AB += (ra[i] - rb[i]) * (ra[i] - rb[i]);

  // if one of the angular number goes below zero,
  // it means it will not have contribution - the same as giving derivation to a constant
  if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0) return 0;

    //Provide recurrence relation
  else if (ax > 0)
    return ((P[0] - ra[0]) *
            electron_repulsive_integral(ra, rb, ax - 1, ay, az, bx, by, bz,
                                        alpha, beta, m + 1) +
            (ax - 1) / 2.0 / alpha *
            electron_repulsive_integral(ra, rb, ax - 2, ay, az, bx, by, bz,
                                        alpha, beta, m) -
            (ax - 1) / 2.0 / alpha * beta / zeta *
            electron_repulsive_integral(ra, rb, ax - 2, ay, az, bx, by, bz,
                                        alpha, beta, m + 1) +
            bx / 2.0 / zeta *
            electron_repulsive_integral(ra, rb, ax - 1, ay, az, bx - 1, by, bz,
                                        alpha, beta, m + 1));

  else if (ay > 0)
    return ((P[1] - ra[1]) *
            electron_repulsive_integral(ra, rb, ax, ay - 1, az, bx, by, bz,
                                        alpha, beta, m + 1) +
            (ay - 1) / 2.0 / alpha *
            electron_repulsive_integral(ra, rb, ax, ay - 2, az, bx, by, bz,
                                        alpha, beta, m) -
            (ay - 1) / 2.0 / alpha * beta / zeta *
            electron_repulsive_integral(ra, rb, ax, ay - 2, az, bx, by, bz,
                                        alpha, beta, m + 1) +
            by / 2.0 / zeta *
            electron_repulsive_integral(ra, rb, ax, ay - 1, az, bx, by - 1, bz,
                                        alpha, beta, m + 1));

  else if (az > 0)
    return ((P[2] - ra[2]) *
            electron_repulsive_integral(ra, rb, ax, ay, az - 1, bx, by, bz,
                                        alpha, beta, m + 1) +
            (az - 1) / 2.0 / alpha *
            electron_repulsive_integral(ra, rb, ax, ay, az - 2, bx, by, bz,
                                        alpha, beta, m) -
            (az - 1) / 2.0 / alpha * beta / zeta *
            electron_repulsive_integral(ra, rb, ax, ay, az - 2, bx, by, bz,
                                        alpha, beta, m + 1) +
            bz / 2.0 / zeta *
            electron_repulsive_integral(ra, rb, ax, ay, az - 1, bx, by, bz - 1,
                                        alpha, beta, m + 1));

  else if (bx > 0)
    return ((P[0] - rb[0]) *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx - 1, by, bz,
                                        alpha, beta, m + 1) +
            (bx - 1) / 2.0 / beta *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx - 2, by, bz,
                                        alpha, beta, m) -
            (bx - 1) / 2.0 / beta * alpha / zeta *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx - 2, by, bz,
                                        alpha, beta, m + 1) +
            ax / 2.0 / zeta *
            electron_repulsive_integral(ra, rb, ax - 1, ay, az, bx - 1, by, bz,
                                        alpha, beta, m + 1));

  else if (by > 0)
    return ((P[1] - rb[1]) *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx, by - 1, bz,
                                        alpha, beta, m + 1) +
            (by - 1) / 2.0 / beta *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx, by - 2, bz,
                                        alpha, beta, m) -
            (by - 1) / 2.0 / beta * alpha / zeta *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx, by - 2, bz,
                                        alpha, beta, m + 1) +
            ay / 2.0 / zeta *
            electron_repulsive_integral(ra, rb, ax, ay - 1, az, bx, by - 1, bz,
                                        alpha, beta, m + 1));

  else if (bz > 0)
    return ((P[2] - rb[2]) *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx, by, bz - 1,
                                        alpha, beta, m + 1) +
            (bz - 1) / 2.0 / beta *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx, by, bz - 2,
                                        alpha, beta, m) -
            (bz - 1) / 2.0 / beta * alpha / zeta *
            electron_repulsive_integral(ra, rb, ax, ay, az, bx, by, bz - 2,
                                        alpha, beta, m + 1) +
            az / 2.0 / zeta *
            electron_repulsive_integral(ra, rb, ax, ay, az - 1, bx, by, bz - 1,
                                        alpha, beta, m + 1));

    //Set the starting point
  else
    return 2.0 * pow(M_PI, 2.5) / alpha / beta / sqrt(zeta) * Boys(xi * AB, m);
}

double
electron_repulsive_integral(const GaussianFunctionPair & pair_1,
                            const GaussianFunctionPair & pair_2) {

  const std::vector<GaussianFunction> AB_pair = expand_function_pair(pair_1);
  const std::vector<GaussianFunction> CD_pair = expand_function_pair(pair_2);

  double result = 0;

  for (const auto & AB_func: AB_pair) {
    for (const auto & CD_func: CD_pair) {
      result += electron_repulsive_integral(AB_func.center.memptr(),
                                            CD_func.center.memptr(),
                                            AB_func.angular[0],
                                            AB_func.angular[1],
                                            AB_func.angular[2],
                                            CD_func.angular[0],
                                            CD_func.angular[1],
                                            CD_func.angular[2],
                                            AB_func.exponent, CD_func.exponent,
                                            0)
                * AB_func.coef * CD_func.coef;

    }
  }

  return result;
}

double
electron_repulsive_integral(const GaussianFunction & A,
                            const GaussianFunction & B,
                            const GaussianFunction & C,
                            const GaussianFunction & D) {
  return electron_repulsive_integral({A, B}, {C, D});
}

double
electron_repulsive_integral(const ERI & eri_info) {
  return electron_repulsive_integral(eri_info.A, eri_info.B, eri_info.C,
                                     eri_info.D);
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
}
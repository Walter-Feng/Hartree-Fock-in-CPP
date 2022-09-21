#include "homebrew.h"

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
    if (abs(PA) < 1e-10 && a - i == 0 && abs(PB) < 1e-19 && b - k + i == 0) {
      result += Binomials(a, i) * Binomials(b, k - i);
    } else if (abs(PA) < 1e-10 && a - i == 0) {
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

// Calculate the transformation coefficient that is needed to complete
// the transformation from two combined gaussian function to the sum
// of a series of gaussian functions
double tranformation_coefficient(int a[3], int b[3], int p[3], double PA[3],
                                 double PB[3], double xi, double AB) {
  int i;
  double result;

  result = 1;

  for (i = 0; i < 3; i++)
    result *= f(p[i], a[i], b[i], PA[i], PB[i]);

  result *= exp(-xi * AB);
  // it might cause confusion, but AB here has already been squared, meaning AB -> |AB|^2

  return result;
}

//Calculate the Overlap Integral of two gaussian function
double SIntegral(double ra[3], double rb[3],
                 int ax, int ay, int az,
                 int bx, int by, int bz,
                 double alpha, double beta) {

  // if one of the angular number goes below zero,
  // it means it will not have contribution - the same as giving derivation to a constant
  if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0) return 0;

    //Provide recurrence relation
  else if (ax > 0)
    return ((alpha * ra[0] + beta * rb[0]) / (alpha + beta) - ra[0]) *
           SIntegral(ra, rb, ax - 1, ay, az, bx, by, bz, alpha, beta) +
           (ax - 1) / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax - 2, ay, az, bx, by, bz, alpha, beta) +
           bx / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha, beta);

  else if (ay > 0)
    return ((alpha * ra[1] + beta * rb[1]) / (alpha + beta) - ra[1]) *
           SIntegral(ra, rb, ax, ay - 1, az, bx, by, bz, alpha, beta) +
           (ay - 1) / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay - 2, az, bx, by, bz, alpha, beta) +
           by / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha, beta);

  else if (az > 0)
    return ((alpha * ra[2] + beta * rb[2]) / (alpha + beta) - ra[2]) *
           SIntegral(ra, rb, ax, ay, az - 1, bx, by, bz, alpha, beta) +
           (az - 1) / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay, az - 2, bx, by, bz, alpha, beta) +
           bz / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha, beta);

  else if (bx > 0)
    return ((alpha * ra[0] + beta * rb[0]) / (alpha + beta) - rb[0]) *
           SIntegral(ra, rb, ax, ay, az, bx - 1, by, bz, alpha, beta) +
           ax / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha, beta) +
           (bx - 1) / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay, az, bx - 2, by, bz, alpha, beta);

  else if (by > 0)
    return ((alpha * ra[1] + beta * rb[1]) / (alpha + beta) - rb[1]) *
           SIntegral(ra, rb, ax, ay, az, bx, by - 1, bz, alpha, beta) +
           ay / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha, beta) +
           (by - 1) / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay, az, bx, by - 2, bz, alpha, beta);

  else if (bz > 0)
    return ((alpha * ra[2] + beta * rb[2]) / (alpha + beta) - rb[2]) *
           SIntegral(ra, rb, ax, ay, az, bx, by, bz - 1, alpha, beta) +
           az / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha, beta) +
           (bz - 1) / 2.0 / (alpha + beta) *
           SIntegral(ra, rb, ax, ay, az, bx, by, bz - 2, alpha, beta);

    //giving the starting point
  else
    return sqrt(M_PI / (alpha + beta)) * M_PI / (alpha + beta) *
           exp(-alpha * beta / (alpha + beta) *
               (pow(ra[0] - rb[0], 2) + pow(ra[1] - rb[1], 2) +
                pow(ra[2] - rb[2], 2)));
}

//Calculate the Coulomb Integral of two gaussian function
double
JIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by,
          int bz, double alpha, double beta, int m) {
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
            JIntegral(ra, rb, ax - 1, ay, az, bx, by, bz, alpha, beta, m + 1) +
            (ax - 1) / 2.0 / alpha *
            JIntegral(ra, rb, ax - 2, ay, az, bx, by, bz, alpha, beta, m) -
            (ax - 1) / 2.0 / alpha * beta / zeta *
            JIntegral(ra, rb, ax - 2, ay, az, bx, by, bz, alpha, beta, m + 1) +
            bx / 2.0 / zeta * JIntegral(ra, rb, ax - 1, ay, az, bx - 1, by, bz,
                                        alpha, beta, m + 1));

  else if (ay > 0)
    return ((P[1] - ra[1]) *
            JIntegral(ra, rb, ax, ay - 1, az, bx, by, bz, alpha, beta, m + 1) +
            (ay - 1) / 2.0 / alpha *
            JIntegral(ra, rb, ax, ay - 2, az, bx, by, bz, alpha, beta, m) -
            (ay - 1) / 2.0 / alpha * beta / zeta *
            JIntegral(ra, rb, ax, ay - 2, az, bx, by, bz, alpha, beta, m + 1) +
            by / 2.0 / zeta * JIntegral(ra, rb, ax, ay - 1, az, bx, by - 1, bz,
                                        alpha, beta, m + 1));

  else if (az > 0)
    return ((P[2] - ra[2]) *
            JIntegral(ra, rb, ax, ay, az - 1, bx, by, bz, alpha, beta, m + 1) +
            (az - 1) / 2.0 / alpha *
            JIntegral(ra, rb, ax, ay, az - 2, bx, by, bz, alpha, beta, m) -
            (az - 1) / 2.0 / alpha * beta / zeta *
            JIntegral(ra, rb, ax, ay, az - 2, bx, by, bz, alpha, beta, m + 1) +
            bz / 2.0 / zeta * JIntegral(ra, rb, ax, ay, az - 1, bx, by, bz - 1,
                                        alpha, beta, m + 1));

  else if (bx > 0)
    return ((P[0] - rb[0]) *
            JIntegral(ra, rb, ax, ay, az, bx - 1, by, bz, alpha, beta, m + 1) +
            (bx - 1) / 2.0 / beta *
            JIntegral(ra, rb, ax, ay, az, bx - 2, by, bz, alpha, beta, m) -
            (bx - 1) / 2.0 / beta * alpha / zeta *
            JIntegral(ra, rb, ax, ay, az, bx - 2, by, bz, alpha, beta, m + 1) +
            ax / 2.0 / zeta * JIntegral(ra, rb, ax - 1, ay, az, bx - 1, by, bz,
                                        alpha, beta, m + 1));

  else if (by > 0)
    return ((P[1] - rb[1]) *
            JIntegral(ra, rb, ax, ay, az, bx, by - 1, bz, alpha, beta, m + 1) +
            (by - 1) / 2.0 / beta *
            JIntegral(ra, rb, ax, ay, az, bx, by - 2, bz, alpha, beta, m) -
            (by - 1) / 2.0 / beta * alpha / zeta *
            JIntegral(ra, rb, ax, ay, az, bx, by - 2, bz, alpha, beta, m + 1) +
            ay / 2.0 / zeta * JIntegral(ra, rb, ax, ay - 1, az, bx, by - 1, bz,
                                        alpha, beta, m + 1));

  else if (bz > 0)
    return ((P[2] - rb[2]) *
            JIntegral(ra, rb, ax, ay, az, bx, by, bz - 1, alpha, beta, m + 1) +
            (bz - 1) / 2.0 / beta *
            JIntegral(ra, rb, ax, ay, az, bx, by, bz - 2, alpha, beta, m) -
            (bz - 1) / 2.0 / beta * alpha / zeta *
            JIntegral(ra, rb, ax, ay, az, bx, by, bz - 2, alpha, beta, m + 1) +
            az / 2.0 / zeta * JIntegral(ra, rb, ax, ay, az - 1, bx, by, bz - 1,
                                        alpha, beta, m + 1));

    //Set the starting point
  else
    return 2.0 * pow(M_PI, 2.5) / alpha / beta / sqrt(zeta) * Boys(xi * AB, m);
}

// The nuclear attraction integral 
double
ZIntegral(double ra[3], double rb[3], double rz[3], int ax, int ay, int az,
          int bx, int by, int bz, double alpha, double beta, int m) {
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
           ZIntegral(ra, rb, rz, ax - 1, ay, az, bx, by, bz, alpha, beta, m) +
           (rz[0] - ra[0] - beta / (alpha + beta) * (rb[0] - ra[0])) *
           ZIntegral(ra, rb, rz, ax - 1, ay, az, bx, by, bz, alpha, beta,
                     m + 1) + (double) (ax - 1) / 2.0 / (alpha + beta) *
                              (ZIntegral(ra, rb, rz, ax - 2, ay, az, bx, by, bz,
                                         alpha, beta, m) -
                               ZIntegral(ra, rb, rz, ax - 2, ay, az, bx, by, bz,
                                         alpha, beta, m + 1)) +
           (double) bx / 2.0 / (alpha + beta) *
           (ZIntegral(ra, rb, rz, ax - 1, ay, az, bx - 1, by, bz, alpha, beta,
                      m) - ZIntegral(ra, rb, rz, ax - 1, ay, az, bx - 1, by, bz,
                                     alpha, beta, m + 1));

  else if (ay > 0)
    return beta / (alpha + beta) * (rb[1] - ra[1]) *
           ZIntegral(ra, rb, rz, ax, ay - 1, az, bx, by, bz, alpha, beta, m) +
           (rz[1] - ra[1] - beta / (alpha + beta) * (rb[1] - ra[1])) *
           ZIntegral(ra, rb, rz, ax, ay - 1, az, bx, by, bz, alpha, beta,
                     m + 1) + (double) (ay - 1) / 2.0 / (alpha + beta) *
                              (ZIntegral(ra, rb, rz, ax, ay - 2, az, bx, by, bz,
                                         alpha, beta, m) -
                               ZIntegral(ra, rb, rz, ax, ay - 2, az, bx, by, bz,
                                         alpha, beta, m + 1)) +
           (double) by / 2.0 / (alpha + beta) *
           (ZIntegral(ra, rb, rz, ax, ay - 1, az, bx, by - 1, bz, alpha, beta,
                      m) - ZIntegral(ra, rb, rz, ax, ay - 1, az, bx, by - 1, bz,
                                     alpha, beta, m + 1));

  else if (az > 0)
    return beta / (alpha + beta) * (rb[2] - ra[2]) *
           ZIntegral(ra, rb, rz, ax, ay, az - 1, bx, by, bz, alpha, beta, m) +
           (rz[2] - ra[2] - beta / (alpha + beta) * (rb[2] - ra[2])) *
           ZIntegral(ra, rb, rz, ax, ay, az - 1, bx, by, bz, alpha, beta,
                     m + 1) + (double) (az - 1) / 2.0 / (alpha + beta) *
                              (ZIntegral(ra, rb, rz, ax, ay, az - 2, bx, by, bz,
                                         alpha, beta, m) -
                               ZIntegral(ra, rb, rz, ax, ay, az - 2, bx, by, bz,
                                         alpha, beta, m + 1)) +
           (double) bz / 2.0 / (alpha + beta) *
           (ZIntegral(ra, rb, rz, ax, ay, az - 1, bx, by, bz - 1, alpha, beta,
                      m) - ZIntegral(ra, rb, rz, ax, ay, az - 1, bx, by, bz - 1,
                                     alpha, beta, m + 1));

  else if (bx > 0)
    return alpha / (alpha + beta) * (ra[0] - rb[0]) *
           ZIntegral(ra, rb, rz, ax, ay, az, bx - 1, by, bz, alpha, beta, m) +
           (rz[0] - rb[0] - alpha / (alpha + beta) * (ra[0] - rb[0])) *
           ZIntegral(ra, rb, rz, ax, ay, az, bx - 1, by, bz, alpha, beta,
                     m + 1) + (double) (bx - 1) / 2.0 / (alpha + beta) *
                              (ZIntegral(ra, rb, rz, ax, ay, az, bx - 2, by, bz,
                                         alpha, beta, m) -
                               ZIntegral(ra, rb, rz, ax, ay, az, bx - 2, by, bz,
                                         alpha, beta, m + 1)) +
           (double) ax / 2.0 / (alpha + beta) *
           (ZIntegral(ra, rb, rz, ax - 1, ay, az, bx - 1, by, bz, alpha, beta,
                      m) - ZIntegral(ra, rb, rz, ax - 1, ay, az, bx - 1, by, bz,
                                     alpha, beta, m + 1));

  else if (by > 0)
    return alpha / (alpha + beta) * (ra[1] - rb[1]) *
           ZIntegral(ra, rb, rz, ax, ay, az, bx, by - 1, bz, alpha, beta, m) +
           (rz[1] - rb[1] - alpha / (alpha + beta) * (ra[1] - rb[1])) *
           ZIntegral(ra, rb, rz, ax, ay, az, bx, by - 1, bz, alpha, beta,
                     m + 1) + (double) (by - 1) / 2.0 / (alpha + beta) *
                              (ZIntegral(ra, rb, rz, ax, ay, az, bx, by - 2, bz,
                                         alpha, beta, m) -
                               ZIntegral(ra, rb, rz, ax, ay, az, bx, by - 2, bz,
                                         alpha, beta, m + 1)) +
           (double) ay / 2.0 / (alpha + beta) *
           (ZIntegral(ra, rb, rz, ax, ay - 1, az, bx, by - 1, bz, alpha, beta,
                      m) - ZIntegral(ra, rb, rz, ax, ay - 1, az, bx, by - 1, bz,
                                     alpha, beta, m + 1));

  else if (bz > 0)
    return alpha / (alpha + beta) * (ra[2] - rb[2]) *
           ZIntegral(ra, rb, rz, ax, ay, az, bx, by, bz - 1, alpha, beta, m) +
           (rz[2] - rb[2] - alpha / (alpha + beta) * (ra[2] - rb[2])) *
           ZIntegral(ra, rb, rz, ax, ay, az, bx, by, bz - 1, alpha, beta,
                     m + 1) + (double) (bz - 1) / 2.0 / (alpha + beta) *
                              (ZIntegral(ra, rb, rz, ax, ay, az, bx, by, bz - 2,
                                         alpha, beta, m) -
                               ZIntegral(ra, rb, rz, ax, ay, az, bx, by, bz - 2,
                                         alpha, beta, m + 1)) +
           (double) az / 2.0 / (alpha + beta) *
           (ZIntegral(ra, rb, rz, ax, ay, az - 1, bx, by, bz - 1, alpha, beta,
                      m) - ZIntegral(ra, rb, rz, ax, ay, az - 1, bx, by, bz - 1,
                                     alpha, beta, m + 1));

  else
    return 2 * M_PI / (alpha + beta) * exp(-alpha * beta / zeta * AB) *
           Boys((alpha + beta) * PC, m);
}
}
#include "integral.h"

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

}


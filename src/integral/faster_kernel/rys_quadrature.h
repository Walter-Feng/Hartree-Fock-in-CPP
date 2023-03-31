#include "integral.h"

namespace hfincpp::integral::faster_kernel::rys_quadrature {

void electron_repulsive_integral_kernel(double * cache,
                                        double C00,
                                        double D00,
                                        double B01,
                                        double B10,
                                        double B00,
                                        int a,
                                        int b,
                                        int stride);

}


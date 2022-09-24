#ifndef HFINCPP_HOMEBREW_INTEGRAL_H
#define HFINCPP_HOMEBREW_INTEGRAL_H

#include "integral.h"

namespace integral::obara_saika {

std::vector<GaussianFunction>
expand_function_pair(const GaussianFunctionPair & pair);

double overlap_integral(const GaussianFunction & A, const GaussianFunction & B);

double overlap_integral(const GaussianFunctionPair & pair);

double electron_repulsive_integral(const GaussianFunction & A,
                                   const GaussianFunction & B,
                                   const GaussianFunction & C,
                                   const GaussianFunction & D);

double electron_repulsive_integral(const GaussianFunctionPair & pair_1,
                                   const GaussianFunctionPair & pair_2);

double electron_repulsive_integral(const ERI & eri_info);
}


#endif //HFINCPP_HOMEBREW_INTEGRAL_H

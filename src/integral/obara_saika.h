#ifndef INTEGRAL_OBARA_SAIKA_H
#define INTEGRAL_OBARA_SAIKA_H

#include "integral.h"

#include "basis/basis.h"

namespace hfincpp::integral::obara_saika {

double overlap_integral(const GaussianFunction & A, const GaussianFunction & B);

double overlap_integral(const GaussianFunctionPair & pair);

arma::mat overlap_integral(const basis::Basis & basis);

arma::mat kinetic_integral(const basis::Basis & basis);

double electron_repulsive_integral(const GaussianFunction & A,
                                   const GaussianFunction & B,
                                   const GaussianFunction & C,
                                   const GaussianFunction & D);

double electron_repulsive_integral(const GaussianFunctionPair & pair_1,
                                   const GaussianFunctionPair & pair_2);

double electron_repulsive_integral(const ERI & eri_info);

arma::mat electron_repulsive_integral(const basis::Basis & basis);

namespace gradient {
arma::cube overlap_integral(const basis::Basis & basis);
}

}


#endif //INTEGRAL_OBARA_SAIKA_H

#ifndef GRADIENT_UTILS_H
#define GRADIENT_UTILS_H

#include <armadillo>

namespace hfincpp::gradient::utils {

arma::cube inverse_r(const arma::mat & xyz);

}

#endif //GRADIENT_UTILS_H

#ifndef GRADIENT_DRIVER_H
#define GRADIENT_DRIVER_H

#include <functional>
#include <armadillo>
#include <json.hpp>

#include "geometry/geometry.h"



namespace hfincpp::gradient {

using EnergyDriver = std::function<double(const geometry::Atoms &)>;
using GradientDriver = std::function<arma::mat(const geometry::Atoms &)>;

template<class Result>
GradientDriver gradient(const Result & result);

nlohmann::json gradient(const nlohmann::json & previous_result);
}

#endif //GRADIENT_DRIVER_H

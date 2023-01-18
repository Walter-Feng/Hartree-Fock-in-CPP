#ifndef GRADIENT_DRIVER_H
#define GRADIENT_DRIVER_H

#include <functional>
#include <armadillo>
#include <json.hpp>

#include "geometry/geometry.h"
#include "basis/basis.h"

namespace hfincpp::gradient {

using EnergyDriver = std::function<double(const geometry::Atoms &)>;
using GradientDriver =
    std::function<std::pair<double, arma::mat>(const geometry::Atoms &)>;

GradientDriver driver(const nlohmann::json & input,
                      const geometry::Atoms & atoms,
                      const basis::Basis & basis,
                      const std::string method);

EnergyDriver energy_driver(const nlohmann::json & input,
                           const geometry::Atoms & atoms,
                           const basis::Basis & basis,
                           const std::string method);
}

#endif //GRADIENT_DRIVER_H

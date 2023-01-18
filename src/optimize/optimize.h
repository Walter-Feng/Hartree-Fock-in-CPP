#ifndef HFINCPP_OPTIMIZE_H
#define HFINCPP_OPTIMIZE_H

#include "geometry/geometry.h"

#include "gradient/driver.h"

namespace hfincpp::optimize {

std::pair<double, geometry::Atoms>
optimize(
    const gradient::EnergyDriver & energy_driver,
    const gradient::GradientDriver & gradient_driver,
    const geometry::Atoms & atoms,
    double initial_step_size = 1e-1,
    double tolerance = 1e-4,
    double gradient_tolerance = 1e-3,
    size_t total_steps = 100,
    int print_level = 1);

nlohmann::json optimize(const nlohmann::json & input,
                        const geometry::Atoms & atoms,
                        const basis::Basis & basis,
                        const std::string & method);

}
#endif //HFINCPP_OPTIMIZE_H

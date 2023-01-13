#ifndef HFINCPP_OPTIMIZE_H
#define HFINCPP_OPTIMIZE_H

#include "geometry/geometry.h"

#include "gradient/driver.h"

namespace hfincpp::optimize {

geometry::Atoms optimize(const gradient::GradientDriver & gradient_driver,
                         const geometry::Atoms & atoms,
                         std::string optimization_method);

}
#endif //HFINCPP_OPTIMIZE_H

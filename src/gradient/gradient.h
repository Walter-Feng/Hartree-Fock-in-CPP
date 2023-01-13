#ifndef HFINCPP_GRADIENT_H
#define HFINCPP_GRADIENT_H

#include <json.hpp>

#include "geometry/geometry.h"
#include "basis/basis.h"

namespace hfincpp::gradient {

nlohmann::json driver(const nlohmann::json & input,
                      const geometry::Atoms & atoms,
                      const basis::Basis & basis,
                      const std::string method);


}

#endif //HFINCPP_GRADIENT_H

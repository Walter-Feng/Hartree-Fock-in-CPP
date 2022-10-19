#ifndef HF_RHF_H
#define HF_RHF_H

#include <json.hpp>
#include "geometry/geometry.h"
#include "basis/basis.h"

namespace hfincpp {
namespace hf {

nlohmann::json rhf(const nlohmann::json & input,
                   const geometry::Atoms & atoms,
                   const basis::Basis & basis);
}
}


#endif //HF_RHF_H

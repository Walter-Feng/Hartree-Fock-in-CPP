#ifndef GEOMETRY_RESOLVE_H
#define GEOMETRY_RESOLVE_H

#include <json.hpp>
#include "geometry.h"
#include "periodic_table.h"

#include "global/error.h"
#include "util/resolve.h"

namespace hfincpp {
namespace util {

template<>
geometry::Atoms resolve(const nlohmann::json & pt);


}
}



#endif //GEOMETRY_RESOLVE_H

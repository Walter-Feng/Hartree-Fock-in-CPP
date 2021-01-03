#ifndef UTIL_RESOLVE_H
#define UTIL_RESOLVE_H

#include <nlohmann/json.hpp>

namespace hfincpp {
namespace util {

template<typename T>
    T resolve(const nlohmann::json & input);

}
}

#endif //UTIL_RESOLVE_H

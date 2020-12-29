#ifndef UTIL_RESOLVE_H
#define UTIL_RESOLVE_H

#include <boost/property_tree/ptree.hpp>
namespace hfincpp {
namespace util {

namespace ptree = boost::property_tree;
template<typename T>
    T resolve(const ptree::ptree & pt);

}
}

#endif //UTIL_RESOLVE_H

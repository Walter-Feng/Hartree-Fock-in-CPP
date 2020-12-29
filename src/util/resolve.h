#ifndef EXE_UTIL_RESOLVE_H
#define EXE_RESOLVE_RESOLVE_H

#include <boost/property_tree/ptree.hpp>
namespace hfincpp {
namespace util {

namespace ptree = boost::property_tree;
template<typename T>
    T resolve(const ptree::ptree & pt);

}
}

#endif //EXE_RESOLVE_RESOLVE_H

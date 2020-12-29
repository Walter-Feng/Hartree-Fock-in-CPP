#ifndef HFINCPP_RUN_H
#define HFINCPP_RUN_H

#include <boost/property_tree/ptree.hpp>

namespace hfincpp {

namespace ptree = boost::property_tree;

ptree::ptree run(const ptree::ptree & input);

}

#endif //HFINCPP_RUN_H

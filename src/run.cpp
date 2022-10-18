#include "run.h"
#include "util/json.h"
#include "geometry/resolve.h"
namespace hfincpp {

nlohmann::json run(const nlohmann::json & input) {

  const geometry::Atoms atoms = util::resolve<geometry::Atoms>(input);

  return {};
}

}

#include "run.h"
#include "util/json.h"
#include "basis/basis.h"
#include "geometry/resolve.h"

#include "hf/rhf.h"

namespace hfincpp {

nlohmann::json run(const nlohmann::json & input) {

  const geometry::Atoms atoms = util::resolve<geometry::Atoms>(input.at("geometry"));
  const std::string basis_string = input.at("basis");
  const basis::Basis basis(atoms, basis_string);
  const std::string method = input.at("method");

  nlohmann::json output;
  output["input"] = input;

  if(method == "rhf") {
    output["output"] = hf::rhf(input, atoms, basis);
  }

  return output;
}

}

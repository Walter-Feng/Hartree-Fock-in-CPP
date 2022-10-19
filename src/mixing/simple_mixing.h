#ifndef MIXING_SIMPLE_MIXING_H
#define MIXING_SIMPLE_MIXING_H

#include "scf/scf.h"

namespace hfincpp::mixing {
template<class State>
struct simple_mixing {
  double alpha;
  State state;

  std::pair<scf::Update<State>, simple_mixing> operator()(const State & state) {
    simple_mixing<State> new_mixer{alpha, state};
    const auto & old_state = state;
    const double alpha_in_function = alpha;
    const scf::Update<State> updater = [old_state, alpha_in_function]
        (const State & state_in_function) -> State {
      const State left = alpha_in_function * state_in_function;
      const State right = (1.0 - alpha_in_function) * old_state;
      return left + right;
    };
    return  {updater, new_mixer};
  }
};
}

#endif //MIXING_SIMPLE_MIXING_H

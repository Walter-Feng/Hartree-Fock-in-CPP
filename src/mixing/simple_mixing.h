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
    const scf::Update<State> updater = [old_state, alpha = this->alpha]
        (const State & state_in_function) -> State {
      const State left = alpha * state_in_function;
      const State right = (1.0 - alpha) * old_state;
      return left + right;
    };
    return  {updater, new_mixer};
  }
};
}

#endif //MIXING_SIMPLE_MIXING_H

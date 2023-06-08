#ifndef N_BODY_SIM_PARTICLE_HPP
#define N_BODY_SIM_PARTICLE_HPP

#pragma once

namespace nbs {
    template <typename Precision, typename Candidate>
    concept Particle = std::floating_point<Precision>
                       && requires (Candidate x) {
                              {
                                  x.position
                                  } -> std::same_as<vec<3, Precision>>;
                          };
}

#endif  // N_BODY_SIM_PARTICLE_HPP

#ifndef N_BODY_SIM_PARTICLE_HPP
#define N_BODY_SIM_PARTICLE_HPP

#pragma once

namespace nbs {
    template <typename Candidate>
    concept Particle = requires (Candidate x) {
                           {
                               x.position
                               } -> std::same_as<vec<3, NBS_PRECISION>&>;
                       };

    template <typename Candidate>
    concept ClusteredParticle
        = Particle<Candidate> && requires (Candidate x) {
                                     {
                                         x.cluster_metadata_idx
                                         } -> std::same_as<size_t&>;
                                 };
}  // namespace nbs

#endif  // N_BODY_SIM_PARTICLE_HPP

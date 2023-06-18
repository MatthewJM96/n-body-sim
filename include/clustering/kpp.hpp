#ifndef N_BODY_SIM_CLUSTERING_KPP_HPP
#define N_BODY_SIM_CLUSTERING_KPP_HPP

#pragma once

#include "particle.hpp"

namespace nbs {
    namespace cluster {
        template <
            size_t                        Dimensions,
            ClusteredParticle<Dimensions> ParticleType,
            KMeansOptions                 Options>
        void
        kpp(const ParticleType* particles,
            IN OUT Cluster<Dimensions, ParticleType>* clusters,
            ui32*                                     seed = nullptr);
    }  // namespace cluster
}  // namespace nbs

#include "kpp.inl"

#endif  // N_BODY_SIM_CLUSTERING_KPP_HPP

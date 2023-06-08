#ifndef N_BODY_SIM_CLUSTERING_KPP_HPP
#define N_BODY_SIM_CLUSTERING_KPP_HPP

#pragma once

#include "particle.hpp"

namespace nbs {
    namespace cluster {
        template <ClusteredParticle ParticleType>
        void
        kpp(ParticleType* particles,
            ui32          particle_count,
            ParticleType* centroids,
            ui32          centroid_count);
    }  // namespace cluster
}  // namespace nbs

#include "kpp.inl"

#endif  // N_BODY_SIM_CLUSTERING_KPP_HPP

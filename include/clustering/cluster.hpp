#ifndef N_BODY_SIM_CLUSTERING_CLUSTER_HPP
#define N_BODY_SIM_CLUSTERING_CLUSTER_HPP

#pragma once

#include "particle.hpp"

namespace nbs {
    namespace cluster {
        template <size_t Dimensions, ClusteredParticle<Dimensions> ParticleType>
        struct Cluster {
            ParticleType centroid;
            size_t       particle_offset;
            size_t       particle_count;
        };
    }  // namespace cluster
}  // namespace nbs

#endif  // N_BODY_SIM_CLUSTERING_CLUSTER_HPP

#ifndef N_BODY_SIM_STATISTICS_AVERAGE_CLUSTER_DISTANCE_HPP
#define N_BODY_SIM_STATISTICS_AVERAGE_CLUSTER_DISTANCE_HPP

#pragma once

#include "clustering/cluster.hpp"
#include "particle.hpp"

namespace nbs {
    namespace statistics {
        template <
            size_t                        Dimensions,
            ClusteredParticle<Dimensions> ParticleType,
            size_t                        ParticleCoun,
            size_t                        ClusterCount>
        f32 calculate_average_cluster_distance(
            ParticleType*                               particles,
            cluster::Cluster<Dimensions, ParticleType>* clusters
        );
    }  // namespace statistics
}  // namespace nbs

#include "average_cluster_distance_inl"

#endif  // N_BODY_SIM_STATISTICS_AVERAGE_CLUSTER_DISTANCE_HPP

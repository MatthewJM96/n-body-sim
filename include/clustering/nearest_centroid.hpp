#ifndef N_BODY_SIM_CLUSTERING_NEAREST_CENTROID_HPP
#define N_BODY_SIM_CLUSTERING_NEAREST_CENTROID_HPP

#pragma once

#include "particle.hpp"

#include "clustering/cluster.hpp"

namespace nbs {
    namespace cluster {
        namespace detail {
            template <Particle ParticleType>
            NearestCentroid nearest_centroid(
                const ParticleType&            particle,
                const ParticleClusterMetadata& particle_metadata,
                const Cluster<ParticleType>*   clusters,
                ui32                           cluster_count
            );
            template <Particle ParticleType>
            NearestCentroid nearest_centroid_from_subset(
                const ParticleType&            particle,
                const ParticleClusterMetadata& particle_metadata,
                const Cluster<ParticleType>*   clusters,
                NearestCentroidList            cluster_subset,
                ui32                           cluster_count
            );

            template <Particle ParticleType>
            std::tuple<NearestCentroid, NearestCentroidList>
            nearest_centroid_and_build_list(
                const ParticleType&            particle,
                const ParticleClusterMetadata& particle_metadata,
                const Cluster<ParticleType>*   clusters,
                ui32                           cluster_count,
                ui32                           subset_count
            );
        };  // namespace detail
    }       // namespace cluster
}  // namespace nbs

#include "nearest_centroid.inl"

#endif  // N_BODY_SIM_CLUSTERING_NEAREST_CENTROID_HPP

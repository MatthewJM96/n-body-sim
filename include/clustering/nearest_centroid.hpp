#ifndef N_BODY_SIM_CLUSTERING_NEAREST_CENTROID_HPP
#define N_BODY_SIM_CLUSTERING_NEAREST_CENTROID_HPP

#pragma once

#include "particle.hpp"

#include "clustering/buffers.hpp"
#include "clustering/cluster.hpp"
#include "clustering/options.hpp"

namespace nbs {
    namespace cluster {
        namespace detail {
            template <ClusteredParticle ParticleType, KMeansOptions Options>
            void nearest_centroid(
                const ParticleType&          particle,
                IN OUT NearestCentroid&      nearest_centroid,
                const Cluster<ParticleType>* clusters
            );
            template <ClusteredParticle ParticleType, KMeansOptions Options>
            void nearest_centroid_from_subset(
                const ParticleType&          particle,
                IN OUT NearestCentroid&      nearest_centroid,
                const Cluster<ParticleType>* clusters,
                NearestCentroidList          cluster_subset
            );

            template <ClusteredParticle ParticleType, KMeansOptions Options>
            void nearest_centroid_and_build_list(
                const ParticleType&          particle,
                OUT NearestCentroid&         nearest_centroid,
                const Cluster<ParticleType>* clusters,
                OUT NearestCentroidList      cluster_subset,
                KMeansBuffers<Options>       buffers
            );
        };  // namespace detail
    }       // namespace cluster
}  // namespace nbs

#include "nearest_centroid.inl"

#endif  // N_BODY_SIM_CLUSTERING_NEAREST_CENTROID_HPP

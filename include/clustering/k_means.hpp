#ifndef N_BODY_SIM_CLUSTERING_K_MEANS_HPP
#define N_BODY_SIM_CLUSTERING_K_MEANS_HPP

#pragma once

#include "particle.hpp"

#include "clustering/buffers.hpp"
#include "clustering/cluster.hpp"
#include "clustering/options.hpp"

namespace nbs {
    namespace cluster {
        template <ClusteredParticle ParticleType, KMeansOptions Options>
        void k_means(
            IN OUT CALLER_DELETE ParticleType* particles,
            IN CALLER_DELETE const Cluster<ParticleType>* initial_clusters,
            OUT CALLER_DELETE Cluster<ParticleType>* final_clusters,
            IN OUT KMeansBuffers<Options> buffers
        );
    }  // namespace cluster
}  // namespace nbs

#include "k_means.inl"

#endif  // N_BODY_SIM_CLUSTERING_K_MEANS_HPP
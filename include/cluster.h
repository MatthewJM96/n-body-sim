#ifndef N_BODY_SIM_CLUSTER_H
#define N_BODY_SIM_CLUSTER_H

#pragma once

#include "particle.hpp"

namespace nbs {
    namespace cluster {
        template <Particle ParticleType>
        struct Cluster {
            ParticleType  centroid;
            ParticleType* particles;
            ui32          particle_count;
        };

        struct KMeansOptions {
            ui32 max_iterations                       = 100;
            ui32 acceptable_changes_per_iteration     = 0;
            bool front_loaded                         = false;
            bool no_approaching_centroid_optimisation = false;
            bool no_centroid_subset_optimisation      = true;

            struct {
                ui32 k_prime = 30;
            } centroid_subset;
        };

        template <Particle ParticleType>
        NBS_PRECISION
        particle_distance_2(const ParticleType& lhs, const ParticleType& rhs);

        template <Particle ParticleType>
        void
        kpp(ParticleType* particles,
            ui32          particle_count,
            ParticleType* centroids,
            ui32          centroid_count);

        template <Particle ParticleType>
        void k_means(
            const Cluster<ParticleType>* initial_clusters,
            ui32                         initial_cluster_count,
            const KMeansOptions&         options,
            OUT Cluster<ParticleType>*& clusters
        );

        namespace impl {
            struct NearestCentroid {
                ui32          idx;
                NBS_PRECISION distance;
            };

            struct ParticleClusterMetadata {
                ui32 initial_cluster_idx;
                ui32 initial_particle_idx;

                NearestCentroid current_cluster;
            };

            struct NearestCentroidList {
                ui32* indices;
            };

            struct NearestCentroidAndList {
                NearestCentroid     centroid;
                NearestCentroidList list;
            };

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
            NearestCentroidAndList nearest_centroid_and_build_list(
                const ParticleType&            particle,
                const ParticleClusterMetadata& particle_metadata,
                const Cluster<ParticleType>*   clusters,
                ui32                           cluster_count,
                ui32                           subset_count
            );
        };  // namespace impl
    }       // namespace cluster
}  // namespace nbs

#include "cluster.inl"

#endif  // N_BODY_SIM_CLUSTER_H

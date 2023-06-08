#include "nearest_centroid.hpp"

template <nbs::ClusteredParticle ParticleType, nbs::cluster::KMeansOptions Options>
void nbs::cluster::k_means(
    IN OUT CALLER_DELETE ParticleType* particles,
    IN CALLER_DELETE const Cluster<ParticleType>* initial_clusters,
    OUT CALLER_DELETE Cluster<ParticleType>* final_clusters,
    IN OUT KMeansBuffers<Options> buffers
) {
    /************
       Set up final clusters for k-means algorithm.
                                    ************/

    // Need to set centroids to initial centroids so as to be useful for the first
    // iteration, in which we have not yet determined our first understanding of where
    // to put the centroids.
    for (size_t cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
        final_clusters[cluster_idx].centroid = initial_clusters[cluster_idx].centroid;
    }

    /************
       Perform k-means algorithm.
                        ************/

    ui32 iterations           = 0;
    ui32 changes_in_iteration = 0;
    do {
        // Complete if max iterations has been reached.
        if (++iterations > Options.max_iterations) break;

        // Each iteration starts at zero changes!
        changes_in_iteration = 0;
        std::fill_n(
            buffers.cluster_modified_in_iteration, Options.cluster_count, false
        );

        //
        // Iterate each particle of the population on which the clusters are being
        // built. For each particle, determine which centroid it is nearest to and add
        // its position to a new centroid which will then take its position as the
        // average position of the associated particles.
        //
        // Iterate each initial cluster and then iterate particles held by that cluster.
        for (ui32 initial_cluster_idx = 0; initial_cluster_idx < Options.cluster_count;
             ++initial_cluster_idx)
        {
            // Get handle on cluster we're looking at.
            const Cluster<ParticleType>& initial_cluster
                = initial_clusters[initial_cluster_idx];

            for (ui32 in_cluster_particle_idx = 0;
                 in_cluster_particle_idx < initial_cluster.particle_count;
                 ++in_cluster_particle_idx)
            {
                ui32 global_particle_idx
                    = initial_cluster.particle_offset + in_cluster_particle_idx;
                detail::NearestCentroid& nearest_centroid
                    = buffers.particle_nearest_centroid[particles[global_particle_idx]
                                                            .cluster_metadata_idx];
                detail::NearestCentroid initial_nearest_centroid = nearest_centroid;

                // If we're front loaded, that means we're starting from no known
                // clusters so just set distance to minimum possible.
                //     This will result in the right behaviour when we search for
                //     the nearest centroid, with calculation performed the first
                //     time over all centroids.
                if constexpr (Options.front_loaded) {
                    nearest_centroid.distance
                        = std::numeric_limits<NBS_PRECISION>::min();
                }

                //
                // Calculate nearest centroid for the current particle, using subset of
                // clusters if we can, and rebuilding that subset if we must.
                //

                if constexpr (Options.centroid_subset_optimisation) {
                    if constexpr (Options.centroid_subset.do_rebuild) {
                        if (iterations == 0)
                            detail::nearest_centroid<ParticleType, Options>(
                                particles
                                    [initial_cluster.particle_offset
                                     + in_cluster_particle_idx],
                                nearest_centroid,
                                final_clusters
                            );
                        else if (iterations == 1) {
                            detail::
                                nearest_centroid_and_build_list<ParticleType, Options>(
                                    particles
                                        [initial_cluster.particle_offset
                                         + in_cluster_particle_idx],
                                    nearest_centroid,
                                    final_clusters,
                                    buffers.nearest_centroids_lists[global_particle_idx]
                                );
                        } else {
                            detail::nearest_centroid_from_subset<ParticleType, Options>(
                                particles
                                    [initial_cluster.particle_offset
                                     + in_cluster_particle_idx],
                                nearest_centroid,
                                final_clusters,
                                buffers.nearest_centroids_lists[global_particle_idx]
                            );
                        }
                    } else {
                        detail::nearest_centroid_from_subset<ParticleType, Options>(
                            particles
                                [initial_cluster.particle_offset
                                 + in_cluster_particle_idx],
                            nearest_centroid,
                            final_clusters,
                            buffers.nearest_centroids_lists[global_particle_idx]
                        );
                    }
                } else {
                    detail::nearest_centroid<ParticleType, Options>(
                        particles
                            [initial_cluster.particle_offset + in_cluster_particle_idx],
                        nearest_centroid,
                        final_clusters
                    );
                }

                // If this is the first particle to join a cluster this round, then set
                // values, otherwise add the new values in.
                if (!buffers.cluster_modified_in_iteration[nearest_centroid.idx]) {
                    final_clusters[nearest_centroid.idx].centroid.position.x
                        = particles
                              [initial_cluster.particle_offset
                               + in_cluster_particle_idx]
                                  .position.x;
                    final_clusters[nearest_centroid.idx].centroid.position.y
                        = particles
                              [initial_cluster.particle_offset
                               + in_cluster_particle_idx]
                                  .position.y;
                    final_clusters[nearest_centroid.idx].centroid.position.z
                        = particles
                              [initial_cluster.particle_offset
                               + in_cluster_particle_idx]
                                  .position.z;

                    final_clusters[nearest_centroid.idx].particle_count = 0;

                    buffers.cluster_modified_in_iteration[nearest_centroid.idx] = true;
                } else {
                    final_clusters[nearest_centroid.idx].centroid.position.x
                        += particles
                               [initial_cluster.particle_offset
                                + in_cluster_particle_idx]
                                   .position.x;
                    final_clusters[nearest_centroid.idx].centroid.position.y
                        += particles
                               [initial_cluster.particle_offset
                                + in_cluster_particle_idx]
                                   .position.y;
                    final_clusters[nearest_centroid.idx].centroid.position.z
                        += particles
                               [initial_cluster.particle_offset
                                + in_cluster_particle_idx]
                                   .position.z;

                    ++(final_clusters[nearest_centroid.idx].particle_count);
                }

                // If the particle has changed cluster, then update changes in
                // iteration and its current cluster index.
                if (initial_nearest_centroid.idx != nearest_centroid.idx) {
                    ++changes_in_iteration;
                }
            }

            // Only look at first cluster in buffer if front loaded.
            if constexpr (Options.front_loaded) break;
        }

        // Using total particles associated with each centroid, calculate the new
        // centroid for that group by taking the average of their positions.
        for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
            final_clusters[cluster_idx].centroid.position
                /= static_cast<NBS_PRECISION>(final_clusters[cluster_idx].particle_count
                );
        }
    } while (changes_in_iteration > Options.acceptable_changes_per_iteration);

    // Once we are done figuring how many particles are in each of the clusters, update
    // the final cluster particle offsets into the underlying particle array.
    size_t curr_offset = 0;
    for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
        final_clusters[cluster_idx].particle_offset = curr_offset;
        curr_offset += final_clusters[cluster_idx].particle_count;
    }

    /************
       Sort particles to their final clusters.
                                     ************/

    auto particle_to_cluster_idx = [&final_clusters, &buffers](const auto& particle) {
        return buffers.particle_nearest_centroid[particle.cluster_metadata_idx].idx;
    };

    std::ranges::sort(
        std::span<ParticleType, Options.particle_count>(
            particles, Options.particle_count
        ),
        std::less<>{},
        particle_to_cluster_idx
    );
}

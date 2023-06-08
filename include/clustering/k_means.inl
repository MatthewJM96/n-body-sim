#include "nearest_centroid.hpp"

template <nbs::Particle ParticleType, nbs::cluster::KMeansOptions Options>
void nbs::cluster::k_means(
    IN CALLER_DELETE const Cluster<ParticleType>* initial_clusters,
    OUT CALLER_DELETE Cluster<ParticleType>* final_clusters,
    IN OUT CALLER_DELETE KMeansBuffers<Options> buffers
) {
    /************
       Set up metadata for k-means algorithm.
                                    ************/

    std::fill_n(
        buffers.particle_metadata,
        Options.particle_count,
        detail::ParticleClusterMetadata{}
    );

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

        ui32 global_particle_idx = 0;
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
                // Set metadata at start of k_means algorithm.
                if (iterations == 0) {
                    buffers.particle_cluster_metadata[global_particle_idx]
                        .initial_particle_idx
                        = in_cluster_particle_idx;
                    buffers.particle_cluster_metadata[global_particle_idx]
                        .initial_cluster_idx
                        = initial_cluster_idx;
                    buffers.particle_cluster_metadata[global_particle_idx]
                        .current_cluster.idx
                        = initial_cluster_idx;

                    // If we're front loaded, that means we're starting from no known
                    // clusters so just set distance to minimum possible.
                    //     This will result in the right behaviour when we search for
                    //     the nearest centroid, with calculation performed the first
                    //     time over all centroids.
                    if constexpr (Options.front_loaded) {
                        buffers.particle_cluster_metadata[global_particle_idx]
                            .current_cluster.distance
                            = std::numeric_limits<NBS_PRECISION>::min();
                    } else {
                        buffers.particle_cluster_metadata[global_particle_idx]
                            .current_cluster.distance
                            = math::distance2(
                                initial_cluster.particles[in_cluster_particle_idx]
                                - initial_cluster.centroid
                            );
                    }
                }

                //
                // Calculate nearest centroid for the current particle, using subset of
                // clusters if we can, and rebuilding that subset if we must.
                //

                detail::NearestCentroid nearest_centroid;
                if constexpr (Options.centroid_subset_optimisation) {
                    if constexpr (Options.centroid_subset.do_rebuild) {
                        if (iterations == 0)
                            nearest_centroid
                                = detail::nearest_centroid<ParticleType, Options>(
                                    initial_cluster.particles[in_cluster_particle_idx],
                                    buffers
                                        .particle_cluster_metadata[global_particle_idx],
                                    final_clusters
                                );
                        else if (iterations == 1) {
                            nearest_centroid = detail::nearest_centroid_and_build_list<
                                ParticleType,
                                Options>(
                                initial_cluster.particles[in_cluster_particle_idx],
                                buffers.particle_cluster_metadata[global_particle_idx],
                                final_clusters,
                                buffers.nearest_centroids_lists[global_particle_idx]
                            );
                        } else {
                            nearest_centroid = detail::nearest_centroid_from_subset<
                                ParticleType,
                                Options>(
                                initial_cluster.particles[in_cluster_particle_idx],
                                buffers.particle_cluster_metadata[global_particle_idx],
                                final_clusters,
                                buffers.nearest_centroids_lists[global_particle_idx]
                            );
                        }
                    } else {
                        nearest_centroid = detail::
                            nearest_centroid_from_subset<ParticleType, Options>(
                                initial_cluster.particles[in_cluster_particle_idx],
                                buffers.particle_cluster_metadata[global_particle_idx],
                                final_clusters,
                                buffers.nearest_centroids_lists[global_particle_idx]
                            );
                    }
                } else {
                    nearest_centroid = detail::nearest_centroid<ParticleType, Options>(
                        initial_cluster.particles[in_cluster_particle_idx],
                        buffers.particle_cluster_metadata[global_particle_idx],
                        final_clusters
                    );
                }

                // If this is the first particle to join a cluster this round, then set
                // values, otherwise add the new values in.
                if (!buffers.cluster_modified_in_iteration[nearest_centroid.idx]) {
                    final_clusters[nearest_centroid.idx].centroid.position.x
                        = initial_cluster.particles[in_cluster_particle_idx].position.x;
                    final_clusters[nearest_centroid.idx].centroid.position.y
                        = initial_cluster.particles[in_cluster_particle_idx].position.y;
                    final_clusters[nearest_centroid.idx].centroid.position.z
                        = initial_cluster.particles[in_cluster_particle_idx].position.z;

                    final_clusters[nearest_centroid.idx].particle_count = 0;

                    buffers.cluster_modified_in_iteration[nearest_centroid.idx] = true;
                } else {
                    final_clusters[nearest_centroid.idx].centroid.position.x
                        += initial_cluster.particles[in_cluster_particle_idx]
                               .position.x;
                    final_clusters[nearest_centroid.idx].centroid.position.y
                        += initial_cluster.particles[in_cluster_particle_idx]
                               .position.y;
                    final_clusters[nearest_centroid.idx].centroid.position.z
                        += initial_cluster.particles[in_cluster_particle_idx]
                               .position.z;

                    ++(final_clusters[nearest_centroid.idx].particle_count);
                }

                // If the particle has changed cluster, then update changes in
                // iteration and its current cluster index.
                if (buffers.particle_cluster_metadata[global_particle_idx]
                        .current_cluster.idx
                    != nearest_centroid.idx)
                {
                    ++changes_in_iteration;
                    buffers.particle_cluster_metadata[global_particle_idx]
                        .current_cluster.idx
                        = nearest_centroid.idx;
                }

                // Update distance to nearest centroid.
                buffers.particle_cluster_metadata[global_particle_idx]
                    .current_cluster.distance
                    = nearest_centroid.distance;

                // Incremement global particle index.
                ++global_particle_idx;
            }

            // Only look at first cluster in buffer if front loaded.
            if (Options.front_loaded) break;
        }

        // Using total particles associated with each centroid, calculate the new
        // centroid for that group by taking the average of their positions.
        for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
            final_clusters[cluster_idx].centroid.position
                /= static_cast<NBS_PRECISION>(final_clusters[cluster_idx].particle_count
                );
        }
    } while (changes_in_iteration > Options.acceptable_changes_per_iteration);

    /************
       Assign particles to their final clusters.
                                     ************/

    // TODO(Matthew): Speed up by minimising copies by considering non-moving and moving
    // particles separately.
    //                    Implementation options include:
    //                        - sorted particle buffers (on distance to centroid, e.g.);
    //                        or,
    //                        - on-demand sorting, taking ownership of previous particle
    //                        buffers and extracting changed particles and placing in
    //                        others.

    // Allocate necessary buffers for clusters.
    for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
        final_clusters[cluster_idx].particles
            = new ParticleType[final_clusters[cluster_idx].particle_count];
    }

    // For each particle, copy it into the appropriate cluster buffer.
    std::array<ui32, Options.cluster_count> particle_cursors = {};
    for (ui32 particle_idx = 0; particle_idx < Options.particle_count; ++particle_idx) {
        ui32 new_cluster_offset
            = particle_cursors[buffers.particle_cluster_metadata[particle_idx]
                                   .current_cluster.idx]++;

        const Cluster<ParticleType>& old_cluster
            = initial_clusters[buffers.particle_cluster_metadata[particle_idx]
                                   .initial_cluster_idx];
        const Cluster<ParticleType>& new_cluster
            = final_clusters[buffers.particle_cluster_metadata[particle_idx]
                                 .current_cluster.idx];

        new_cluster.particles[new_cluster_offset]
            = old_cluster.particles[buffers.particle_cluster_metadata[particle_idx]
                                        .initial_particle_idx];
    }
}

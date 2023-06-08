#include "nearest_centroid.hpp"

template <nbs::Particle ParticleType, nbs::cluster::KMeansOptions Options>
void nbs::cluster::k_means(
    IN CALLER_DELETE const Cluster<ParticleType>* initial_clusters,
    OUT CALLER_DELETE Cluster<ParticleType>* clusters,
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
    std::fill_n(buffers.cluster_modified_in_iteration, Options.cluster_count, false);
    if constexpr (Options.centroid_subset_optimisation) {
        std::fill_n(
            buffers.nearest_centroid_lists,
            Options.particle_count,
            detail::NearestCentroidList{}
        );
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
            const Cluster<ParticleType>& cluster = clusters[initial_cluster_idx];

            for (ui32 in_cluster_particle_idx = 0;
                 in_cluster_particle_idx
                 < initial_clusters[initial_cluster_idx].particle_count;
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
                    // Additionally, if we want to avoid the approaching centroid
                    // optimisation, then we should likewise set the distance to the
                    // minimum possible to force distance calculations for all
                    // centroids.
                    if (Options.front_loaded
                        || !Options.approaching_centroid_optimisation)
                    {
                        buffers.particle_cluster_metadata[global_particle_idx]
                            .current_cluster.distance
                            = std::numeric_limits<NBS_PRECISION>::min();
                    } else {
                        buffers.particle_cluster_metadata[global_particle_idx]
                            .current_cluster.distance
                            = math::distance2(
                                cluster.particles[in_cluster_particle_idx]
                                - cluster.centroid
                            );
                    }
                }

                detail::NearestCentroid nearest_centroid;
                if constexpr (Options.centroid_subset_optimisation) {
                    if (iterations == 0)
                        nearest_centroid
                            = detail::nearest_centroid<ParticleType, Options>(
                                cluster.particles[in_cluster_particle_idx],
                                buffers.particle_cluster_metadata[global_particle_idx],
                                clusters
                            );
                    else if (iterations == 1) {
                        nearest_centroid = detail::
                            nearest_centroid_and_build_list<ParticleType, Options>(
                                cluster.particles[in_cluster_particle_idx],
                                buffers.particle_cluster_metadata[global_particle_idx],
                                clusters,
                                buffers.nearest_centroids_lists[global_particle_idx]
                            );
                    } else {
                        nearest_centroid = detail::
                            nearest_centroid_from_subset<ParticleType, Options>(
                                cluster.particles[in_cluster_particle_idx],
                                buffers.particle_cluster_metadata[global_particle_idx],
                                clusters,
                                buffers.nearest_centroids_lists[global_particle_idx]
                            );
                    }
                } else {
                    nearest_centroid = detail::nearest_centroid<ParticleType, Options>(
                        cluster.particles[in_cluster_particle_idx],
                        buffers.particle_cluster_metadata[global_particle_idx],
                        clusters
                    );
                }

                // If this is the first particle to join a cluster this round, then set
                // values, otherwise add the new values in.
                if (!buffers.cluster_modified_in_iteration[initial_cluster_idx]) {
                    clusters[nearest_centroid.idx].centroid.x
                        = cluster.particles[in_cluster_particle_idx].x;
                    clusters[nearest_centroid.idx].centroid.y
                        = cluster.particles[in_cluster_particle_idx].y;
                    clusters[nearest_centroid.idx].centroid.z
                        = cluster.particles[in_cluster_particle_idx].z;

                    clusters[nearest_centroid.idx].particle_count = 0;

                    buffers.cluster_modified_in_iteration[initial_cluster_idx] = true;
                } else {
                    clusters[nearest_centroid.idx].centroid.x
                        += cluster.particles[in_cluster_particle_idx].x;
                    clusters[nearest_centroid.idx].centroid.y
                        += cluster.particles[in_cluster_particle_idx].y;
                    clusters[nearest_centroid.idx].centroid.z
                        += cluster.particles[in_cluster_particle_idx].z;

                    ++(clusters[nearest_centroid.idx].particle_count);
                }

                // If the particle has changed cluster particleship, then update changes
                // in iteration and its current cluster index.
                if (buffers.particle_cluster_metadata[global_particle_idx]
                        .current_cluster.idx
                    != nearest_centroid.idx)
                {
                    ++changes_in_iteration;
                    buffers.particle_cluster_metadata[global_particle_idx]
                        .current_cluster.idx
                        = nearest_centroid.idx;
                }

                // No matter what, the distance to the centroid has very likely changed,
                // so we should update it! However, if we want to avoid the approaching
                // centroid optimisation, then we should set the distance to the minimum
                // possible to force distance calculations for all centroids.
                if (!Options.approaching_centroid_optimisation) {
                    buffers.particle_cluster_metadata[global_particle_idx]
                        .current_cluster.distance
                        = std::numeric_limits<NBS_PRECISION>::min();
                } else {
                    buffers.particle_cluster_metadata[global_particle_idx]
                        .current_cluster.distance
                        = nearest_centroid.distance;
                }

                // Incremement global particle index.
                ++global_particle_idx;
            }

            // Only look at first cluster in buffer if front loaded.
            if (Options.front_loaded) break;
        }

        // Using total particles associated with each centroid, calculate the new
        // centroid for that group by taking the average of their positions.
        for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
            clusters[cluster_idx].centroid.x /= clusters[cluster_idx].particle_count;
            clusters[cluster_idx].centroid.y /= clusters[cluster_idx].particle_count;
            clusters[cluster_idx].centroid.z /= clusters[cluster_idx].particle_count;
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
        clusters[cluster_idx].particles
            = new ParticleType[clusters[cluster_idx].particle_count];
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
            = clusters[buffers.particle_cluster_metadata[particle_idx]
                           .current_cluster.idx];

        new_cluster.particles[new_cluster_offset]
            = old_cluster.particles[buffers.particle_cluster_metadata[particle_idx]
                                        .initial_particle_idx];
    }
}

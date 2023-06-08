template <nbs::Particle ParticleType>
NBS_PRECISION
nbs::cluster::particle_distance_2(const ParticleType& lhs, const ParticleType& rhs) {
    return (lhs.x - rhs.x) * (lhs.x - rhs.x) + (lhs.y - rhs.y) * (lhs.y - rhs.y)
           + (lhs.z - rhs.z) * (lhs.z - rhs.z);
}

template <nbs::Particle ParticleType>
void nbs::cluster::kpp(
    ParticleType* particles,
    ui32          particle_count,
    ParticleType* centroids,
    ui32          centroid_count
) {
    /************
       Set up metadata for k++ algorithm.
                                ************/

    ui32  number_of_particles_chosen  = 1;
    bool* particle_chosen_as_centroid = new bool[particle_count];

    NBS_PRECISION* cumulative_distance_2s = new NBS_PRECISION[particle_count];
    NBS_PRECISION  total_distance_2       = 0.0;

    std::default_random_engine generator;

    /************
       Make initial choice of a centroid.
                                ************/

    {
        std::uniform_int_distribution<ui32> distribution(0, particle_count - 1);
        ui32                                initial_choice = distribution(generator);

        centroids[0]                                = particles[initial_choice];
        particle_chosen_as_centroid[initial_choice] = true;
    }

    /************
       Perform k++ algorithm for each subsequent centroid.
                                                 ************/

    for (ui32 centroid_idx = 1; centroid_idx < centroid_count; ++centroid_idx) {
        //
        // For each particle of the dataset not so far chosen as a centroid, determine
        // the minimum distance to a chosen centroid and select one of those particles
        // to be the next centroid with probability proportional to distance^2 from
        // nearest centroid.
        //

        // Calculate distances for each particle not so far chosen as a centroid.
        for (ui32 particle_idx = 0; particle_idx < particle_count; ++particle_idx) {
            // Skip particles already chosen as a centroid.
            if (particle_chosen_as_centroid[particle_idx]) {
                cumulative_distance_2s[particle_idx]
                    = cumulative_distance_2s[particle_idx - 1];
                continue;
            }

            // Calculate distance to nearest chosen centroid.
            NBS_PRECISION minimum_distance_2_to_chosen_centroids
                = std::numeric_limits<NBS_PRECISION>::max();
            for (ui32 chosen_centroid_idx = 0; chosen_centroid_idx < centroid_idx;
                 ++chosen_centroid_idx)
            {
                NBS_PRECISION distance_2_to_chosen_centroid = particle_distance_2(
                    particles[particle_idx], centroids[chosen_centroid_idx]
                );

                if (distance_2_to_chosen_centroid
                    < minimum_distance_2_to_chosen_centroids)
                {
                    minimum_distance_2_to_chosen_centroids
                        = distance_2_to_chosen_centroid;
                }
            }

            // Place the calculated distance into metadata.
            if (particle_idx == 0) {
                cumulative_distance_2s[particle_idx]
                    = minimum_distance_2_to_chosen_centroids;
            } else {
                cumulative_distance_2s[particle_idx]
                    = cumulative_distance_2s[particle_idx - 1]
                      + minimum_distance_2_to_chosen_centroids;
            }
            total_distance_2 += minimum_distance_2_to_chosen_centroids;
        }

        // Select a particle to be next chosen centroid with probability proportional to
        // distance^2 to nearest existing centroid.
        std::uniform_real_distribution<NBS_PRECISION> distribution(
            0.0, total_distance_2
        );
        NBS_PRECISION choice = distribution(generator);

        for (ui32 particle_idx = 0; particle_idx < particle_count; ++particle_idx) {
            if (cumulative_distance_2s[particle_idx] > choice) {
                centroids[centroid_idx]                   = particles[particle_idx];
                particle_chosen_as_centroid[particle_idx] = true;
            }
        }
    }

    /************
       Clean-up.
       ************/

    delete[] particle_chosen_as_centroid;
    delete[] cumulative_distance_2s;
}

template <nbs::Particle ParticleType>
nbs::cluster::impl::NearestCentroid nbs::cluster::impl::nearest_centroid(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters,
    ui32                           cluster_count
) {
    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    NBS_PRECISION new_distance_2_to_current_cluster = particle_distance_2(
        particle, clusters[particle_metadata.current_cluster.idx].centroid
    );
    if (new_distance_2_to_current_cluster < particle_metadata.current_cluster.distance)
    {
        return NearestCentroid{ particle_metadata.current_cluster.idx,
                                new_distance_2_to_current_cluster };
    }

    NearestCentroid nearest_centroid
        = { particle_metadata.current_cluster.idx, new_distance_2_to_current_cluster };

    for (ui32 cluster_idx = 0; cluster_idx < cluster_count; ++cluster_idx) {
        if (cluster_idx == particle_metadata.current_cluster.idx) continue;

        NBS_PRECISION centroid_distance_2
            = particle_distance_2(particle, clusters[cluster_idx].centroid);

        if (centroid_distance_2 < nearest_centroid.distance) {
            nearest_centroid.idx      = cluster_idx;
            nearest_centroid.distance = centroid_distance_2;
        }
    }

    return nearest_centroid;
}

template <nbs::Particle ParticleType>
nbs::cluster::impl::NearestCentroid nbs::cluster::impl::nearest_centroid_from_subset(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters,
    impl::NearestCentroidList      cluster_subset,
    ui32                           cluster_count
) {
    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    NBS_PRECISION new_distance_2_to_current_cluster = particle_distance_2(
        particle, clusters[particle_metadata.current_cluster.idx].centroid
    );
    if (new_distance_2_to_current_cluster < particle_metadata.current_cluster.distance)
    {
        return NearestCentroid{ particle_metadata.current_cluster.idx,
                                new_distance_2_to_current_cluster };
    }

    NearestCentroid nearest_centroid
        = { particle_metadata.current_cluster.idx, new_distance_2_to_current_cluster };

    for (ui32 cluster_subset_idx = 0; cluster_subset_idx < cluster_count;
         ++cluster_subset_idx)
    {
        ui32 cluster_idx = cluster_subset.indices[cluster_subset_idx];

        if (cluster_idx == particle_metadata.current_cluster.idx) continue;

        NBS_PRECISION centroid_distance_2
            = particle_distance_2(particle, clusters[cluster_idx].centroid);

        if (centroid_distance_2 < nearest_centroid.distance) {
            nearest_centroid.idx      = cluster_idx;
            nearest_centroid.distance = centroid_distance_2;
        }
    }

    return nearest_centroid;
}

#include <map>

template <nbs::Particle ParticleType>
nbs::cluster::impl::NearestCentroidAndList
nbs::cluster::impl::nearest_centroid_and_build_list(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters,
    ui32                           cluster_count,
    ui32                           subset_count
) {
    // Optimisation by making a subset of centroids to consider for a given particle.
    //     This is based on the paper "Faster k-means Cluster Estimation" by Khandelwal
    //     S., Awekar A.
    std::multimap<NBS_PRECISION, ui32> centroid_list;

    NearestCentroidAndList nearest_centroid_and_list = {
        NearestCentroid{
                        particle_metadata.current_cluster.idx,
                        particle_distance_2(
                particle, clusters[particle_metadata.current_cluster.idx].centroid
            ) },
        NearestCentroidList{ new ui32[subset_count] }
    };

    for (ui32 cluster_idx = 0; cluster_idx < cluster_count; ++cluster_idx) {
        if (cluster_idx == particle_metadata.current_cluster.idx) continue;

        NBS_PRECISION centroid_distance_2
            = particle_distance_2(particle, clusters[cluster_idx].centroid);

        centroid_list.insert(std::make_pair(centroid_distance_2, cluster_idx));

        if (centroid_distance_2 < nearest_centroid_and_list.centroid.distance) {
            nearest_centroid_and_list.centroid.idx      = cluster_idx;
            nearest_centroid_and_list.centroid.distance = centroid_distance_2;
        }
    }

    auto it = centroid_list.begin();
    for (ui32 idx = 0; idx < subset_count; ++idx) {
        nearest_centroid_and_list.list.indices[idx] = *it++;
    }

    return nearest_centroid_and_list;
}

#include <cstring>

template <nbs::Particle ParticleType>
void nbs::cluster::k_means(
    const Cluster<ParticleType>* initial_clusters,
    ui32                         cluster_count,
    ui32                         particle_count,
    const KMeansOptions&         options,
    OUT Cluster<ParticleType>*& clusters
) {
    /************
       Set up buffer of new clusters produced.
                                     ************/

    clusters = new Cluster<ParticleType>[cluster_count];
    std::memcpy(
        clusters, initial_clusters, sizeof(Cluster<ParticleType>) * cluster_count
    );

    /************
       Set up metadata for k-means algorithm.
                                    ************/

    impl::ParticleClusterMetadata* particle_cluster_metadata
        = new impl::ParticleClusterMetadata[particle_count];
    // TODO(Matthew): Do we template parameterise the centroid_subset_optimisation flag
    // to reduce memory usage if it is disabled?
    [[maybe_unused]] impl::NearestCentroidList* nearest_centroids_lists
        = new impl::NearestCentroidList[particle_count];
    bool* cluster_modified_in_iteration = new bool[cluster_count](false);

    /************
       Perform k-means algorithm.
                        ************/

    ui32 iterations           = 0;
    ui32 changes_in_iteration = 0;
    do {
        // Complete if max iterations has been reached.
        if (++iterations > options.max_iterations) break;

        // Each iteration starts at zero changes!
        changes_in_iteration = 0;
        std::fill_n(cluster_modified_in_iteration, cluster_count, false);

        //
        // Iterate each particle of the population on which the clusters are being
        // built. For each particle, determine which centroid it is nearest to and add
        // its position to a new centroid which will then take its position as the
        // average position of the associated particles.
        //

        ui32 global_particle_idx = 0;
        // Iterate each initial cluster and then iterate particles held by that cluster.
        for (ui32 initial_cluster_idx = 0; initial_cluster_idx < cluster_count;
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
                    particle_cluster_metadata[global_particle_idx].initial_particle_idx
                        = in_cluster_particle_idx;
                    particle_cluster_metadata[global_particle_idx].initial_cluster_idx
                        = initial_cluster_idx;
                    particle_cluster_metadata[global_particle_idx].current_cluster.idx
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
                    if (options.front_loaded
                        || options.no_approaching_centroid_optimisation)
                    {
                        particle_cluster_metadata[global_particle_idx]
                            .current_cluster.distance
                            = std::numeric_limits<NBS_PRECISION>::min();
                    } else {
                        particle_cluster_metadata[global_particle_idx]
                            .current_cluster.distance
                            = particle_distance_2(
                                cluster.particles[in_cluster_particle_idx],
                                cluster.centroid
                            );
                    }
                }

                impl::NearestCentroid nearest_centroid;
                if (options.no_centroid_subset_optimisation || iterations == 0) {
                    nearest_centroid = impl::nearest_centroid(
                        cluster.particles[in_cluster_particle_idx],
                        particle_cluster_metadata[global_particle_idx],
                        clusters,
                        cluster_count
                    );
                } else {
                    if (iterations == 1) {
                        impl::NearestCentroidAndList nearest_centroid_and_list
                            = impl::nearest_centroid_and_build_list(
                                cluster.particles[in_cluster_particle_idx],
                                particle_cluster_metadata[global_particle_idx],
                                clusters,
                                cluster_count,
                                options.centroid_subset.k_prime
                            );

                        nearest_centroid = nearest_centroid_and_list.centroid;
                        nearest_centroids_lists[global_particle_idx]
                            = nearest_centroid_and_list.list;
                    } else {
                        nearest_centroid = impl::nearest_centroid_from_subset(
                            cluster.particles[in_cluster_particle_idx],
                            particle_cluster_metadata[global_particle_idx],
                            clusters,
                            nearest_centroids_lists[global_particle_idx],
                            cluster_count
                        );
                    }
                }

                // If this is the first particle to join a cluster this round, then set
                // values, otherwise add the new values in.
                if (!cluster_modified_in_iteration[initial_cluster_idx]) {
                    clusters[nearest_centroid.idx].centroid.x
                        = cluster.particles[in_cluster_particle_idx].x;
                    clusters[nearest_centroid.idx].centroid.y
                        = cluster.particles[in_cluster_particle_idx].y;
                    clusters[nearest_centroid.idx].centroid.z
                        = cluster.particles[in_cluster_particle_idx].z;

                    clusters[nearest_centroid.idx].particle_count = 0;

                    cluster_modified_in_iteration[initial_cluster_idx] = true;
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
                if (particle_cluster_metadata[global_particle_idx].current_cluster.idx
                    != nearest_centroid.idx)
                {
                    ++changes_in_iteration;
                    particle_cluster_metadata[global_particle_idx].current_cluster.idx
                        = nearest_centroid.idx;
                }

                // No matter what, the distance to the centroid has very likely changed,
                // so we should update it! However, if we want to avoid the approaching
                // centroid optimisation, then we should set the distance to the minimum
                // possible to force distance calculations for all centroids.
                if (options.no_approaching_centroid_optimisation) {
                    particle_cluster_metadata[global_particle_idx]
                        .current_cluster.distance
                        = std::numeric_limits<NBS_PRECISION>::min();
                } else {
                    particle_cluster_metadata[global_particle_idx]
                        .current_cluster.distance
                        = nearest_centroid.distance;
                }

                // Incremement global particle index.
                ++global_particle_idx;
            }

            // Only look at first cluster in buffer if front loaded.
            if (options.front_loaded) break;
        }

        // Using total particles associated with each centroid, calculate the new
        // centroid for that group by taking the average of their positions.
        for (ui32 cluster_idx = 0; cluster_idx < cluster_count; ++cluster_idx) {
            clusters[cluster_idx].centroid.x /= clusters[cluster_idx].particle_count;
            clusters[cluster_idx].centroid.y /= clusters[cluster_idx].particle_count;
            clusters[cluster_idx].centroid.z /= clusters[cluster_idx].particle_count;
        }
    } while (changes_in_iteration > options.acceptable_changes_per_iteration);

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
    for (ui32 cluster_idx = 0; cluster_idx < cluster_count; ++cluster_idx) {
        clusters[cluster_idx].particles
            = new ParticleType[clusters[cluster_idx].particle_count];
    }

    // For each particle, copy it into the appropriate cluster buffer.
    ui32* particle_cursors = new ui32[cluster_count](0);
    for (ui32 particle_idx = 0; particle_idx < particle_count; ++particle_idx) {
        ui32 new_cluster_offset
            = particle_cursors[particle_cluster_metadata[particle_idx]
                                   .current_cluster.idx]++;

        const Cluster<ParticleType>& old_cluster
            = initial_clusters[particle_cluster_metadata[particle_idx]
                                   .initial_cluster_idx];
        const Cluster<ParticleType>& new_cluster
            = clusters[particle_cluster_metadata[particle_idx].current_cluster.idx];

        new_cluster.particles[new_cluster_offset]
            = old_cluster.particles[particle_cluster_metadata[particle_idx]
                                        .initial_particle_idx];
    }

    /************
       Clean-up.
       ************/

    delete[] particle_cluster_metadata;
    delete[] cluster_modified_in_iteration;
}

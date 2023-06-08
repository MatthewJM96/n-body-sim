template <
    size_t                             Dimensions,
    nbs::ClusteredParticle<Dimensions> ParticleType,
    nbs::cluster::KMeansOptions        Options>
void nbs::cluster::kpp(
    const ParticleType* particles, IN OUT Cluster<Dimensions, ParticleType>* clusters
) {
    /************
       Set up metadata for k++ algorithm.
                                ************/

    bool* particle_chosen_as_centroid = new bool[Options.particle_count]{};

    NBS_PRECISION* cumulative_distance_2s = new NBS_PRECISION[Options.particle_count]{};

    std::default_random_engine generator;

    /************
       Make initial choice of a centroid.
                                ************/

    {
        std::uniform_int_distribution<ui32> distribution(0, Options.particle_count - 1);
        ui32                                initial_choice = distribution(generator);

        clusters[0].centroid                        = particles[initial_choice];
        particle_chosen_as_centroid[initial_choice] = true;
    }

    /************
       Perform k++ algorithm for each subsequent centroid.
                                                 ************/

    for (ui32 cluster_idx = 1; cluster_idx < Options.cluster_count; ++cluster_idx) {
        NBS_PRECISION total_distance_2 = 0.0;

        //
        // For each particle of the dataset not so far chosen as a centroid, determine
        // the minimum distance to a chosen centroid and select one of those particles
        // to be the next centroid with probability proportional to distance^2 from
        // nearest centroid.
        //

        // Calculate distances for each particle not so far chosen as a centroid.
        for (ui32 particle_idx = 0; particle_idx < Options.particle_count;
             ++particle_idx)
        {
            // Skip particles already chosen as a centroid.
            if (particle_chosen_as_centroid[particle_idx]) {
                if (particle_idx == 0) {
                    cumulative_distance_2s[0] = 0.0;
                } else {
                    cumulative_distance_2s[particle_idx]
                        = cumulative_distance_2s[particle_idx - 1];
                }
                continue;
            }

            // Calculate distance to nearest chosen centroid.
            NBS_PRECISION minimum_distance_2_to_chosen_centroids
                = std::numeric_limits<NBS_PRECISION>::max();
            for (ui32 chosen_cluster_idx = 0; chosen_cluster_idx < cluster_idx;
                 ++chosen_cluster_idx)
            {
                NBS_PRECISION distance_2_to_chosen_centroid = math::distance2(
                    particles[particle_idx].position,
                    clusters[chosen_cluster_idx].centroid.position
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

        for (ui32 particle_idx = 0; particle_idx < Options.particle_count;
             ++particle_idx)
        {
            if (cumulative_distance_2s[particle_idx] > choice) {
                clusters[cluster_idx].centroid            = particles[particle_idx];
                particle_chosen_as_centroid[particle_idx] = true;
                break;
            }
        }
    }

    /************
       Clean-up.
       ************/

    delete[] particle_chosen_as_centroid;
    delete[] cumulative_distance_2s;
}

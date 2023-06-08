template <nbs::ClusteredParticle ParticleType>
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
                NBS_PRECISION distance_2_to_chosen_centroid = math::distance2(
                    particles[particle_idx] - centroids[chosen_centroid_idx]
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

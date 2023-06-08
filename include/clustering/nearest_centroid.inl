template <
    size_t                             Dimensions,
    nbs::ClusteredParticle<Dimensions> ParticleType,
    nbs::cluster::KMeansOptions        Options>
void nbs::cluster::detail::nearest_centroid(
    const ParticleType&                      particle,
    IN OUT NearestCentroid&                  nearest_centroid,
    const Cluster<Dimensions, ParticleType>* clusters
) {
    NBS_PRECISION new_distance_2_to_current_cluster = math::distance2(
        particle.position, clusters[nearest_centroid.idx].centroid.position
    );

    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    if constexpr (Options.approaching_centroid_optimisation) {
        if (new_distance_2_to_current_cluster < nearest_centroid.distance) {
            nearest_centroid.distance = new_distance_2_to_current_cluster;

            return;
        }
    }

    nearest_centroid.distance = new_distance_2_to_current_cluster;

    // For each centroid, consider if it is closer than the current centroid.
    for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
        if (cluster_idx == nearest_centroid.idx) continue;

        NBS_PRECISION centroid_distance_2 = math::distance2(
            particle.position, clusters[cluster_idx].centroid.position
        );

        if (centroid_distance_2 < nearest_centroid.distance) {
            nearest_centroid.idx      = cluster_idx;
            nearest_centroid.distance = centroid_distance_2;
        }
    }
}

template <
    size_t                             Dimensions,
    nbs::ClusteredParticle<Dimensions> ParticleType,
    nbs::cluster::KMeansOptions        Options>
void nbs::cluster::detail::nearest_centroid_from_subset(
    const ParticleType&                      particle,
    IN OUT NearestCentroid&                  nearest_centroid,
    const Cluster<Dimensions, ParticleType>* clusters,
    detail::NearestCentroidList              cluster_subset
) {
    NBS_PRECISION new_distance_2_to_current_cluster = math::distance2(
        particle.position, clusters[nearest_centroid.idx].centroid.position
    );

    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    if constexpr (Options.approaching_centroid_optimisation) {
        if (new_distance_2_to_current_cluster < nearest_centroid.distance) {
            nearest_centroid.distance = new_distance_2_to_current_cluster;

            return;
        }
    }

    nearest_centroid.distance = new_distance_2_to_current_cluster;

    // For each centroid, consider if it is closer than the current centroid.
    for (ui32 cluster_subset_idx = 0;
         cluster_subset_idx < Options.centroid_subset.k_prime;
         ++cluster_subset_idx)
    {
        ui32 cluster_idx = cluster_subset.indices[cluster_subset_idx];

        if (cluster_idx == nearest_centroid.idx) continue;

        NBS_PRECISION centroid_distance_2 = math::distance2(
            particle.position, clusters[cluster_idx].centroid.position
        );

        if (centroid_distance_2 < nearest_centroid.distance) {
            nearest_centroid.idx      = cluster_idx;
            nearest_centroid.distance = centroid_distance_2;
        }
    }
}

template <
    size_t                             Dimensions,
    nbs::ClusteredParticle<Dimensions> ParticleType,
    nbs::cluster::KMeansOptions        Options>
void nbs::cluster::detail::nearest_centroid_and_build_list(
    const ParticleType&                      particle,
    OUT NearestCentroid&                     nearest_centroid,
    const Cluster<Dimensions, ParticleType>* clusters,
    OUT detail::NearestCentroidList& cluster_subset,
    KMeansBuffers<Options>           buffers
) {
    // Optimisation by making a subset of centroids to consider for a given particle.
    //     This is based on the paper "Faster k-means Cluster Estimation" by Khandelwal
    //     S., Awekar A.

    // Place indices in order.
    for (ui32 i = 0; i < Options.cluster_count; ++i)
        buffers.nearest_centroid_indices[i] = i;

    // Transformer from index to distance to particle building the cluster subset for.
    auto index_to_distance = [&particle, &clusters](ui32 idx) {
        return math::distance2(
            particle.position.position, clusters[idx].centroid.position.position
        );
    };

    // Sort indices according to distance to particle.
    std::ranges::sort(
        buffers.nearest_centroid_indices, std::less<>{}, index_to_distance
    );

    // Set cluster subset.
    for (ui32 idx = 0; idx < Options.centroid_subset.k_prime; ++idx) {
        cluster_subset.indices[idx] = buffers.nearest_centroid_indices[idx];
    }

    nearest_centroid.idx      = buffers.nearest_centroid_indices[0];
    nearest_centroid.distance = math::distance2(
        particle, clusters[buffers.nearest_centroid_indices[0]].centroid.position
    );
}

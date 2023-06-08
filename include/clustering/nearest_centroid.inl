template <nbs::Particle ParticleType, nbs::cluster::KMeansOptions Options>
nbs::cluster::detail::NearestCentroid nbs::cluster::detail::nearest_centroid(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters
) {
    NBS_PRECISION new_distance_2_to_current_cluster = math::distance2(
        particle - clusters[particle_metadata.current_cluster.idx].centroid
    );

    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    if constexpr (Options.approaching_centroid_optimisation) {
        if (new_distance_2_to_current_cluster
            < particle_metadata.current_cluster.distance)
        {
            return NearestCentroid{ particle_metadata.current_cluster.idx,
                                    new_distance_2_to_current_cluster };
        }
    }

    NearestCentroid nearest_centroid
        = { particle_metadata.current_cluster.idx, new_distance_2_to_current_cluster };

    // For each centroid, consider if it is closer than the current centroid.
    for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
        if (cluster_idx == particle_metadata.current_cluster.idx) continue;

        NBS_PRECISION centroid_distance_2
            = math::distance2(particle - clusters[cluster_idx].centroid);

        if (centroid_distance_2 < nearest_centroid.distance) {
            nearest_centroid.idx      = cluster_idx;
            nearest_centroid.distance = centroid_distance_2;
        }
    }

    return nearest_centroid;
}

template <nbs::Particle ParticleType, nbs::cluster::KMeansOptions Options>
nbs::cluster::detail::NearestCentroid
nbs::cluster::detail::nearest_centroid_from_subset(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters,
    detail::NearestCentroidList    cluster_subset
) {
    NBS_PRECISION new_distance_2_to_current_cluster = math::distance2(
        particle - clusters[particle_metadata.current_cluster.idx].centroid
    );

    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    if constexpr (Options.approaching_centroid_optimisation) {
        if (new_distance_2_to_current_cluster
            < particle_metadata.current_cluster.distance)
        {
            return NearestCentroid{ particle_metadata.current_cluster.idx,
                                    new_distance_2_to_current_cluster };
        }
    }

    NearestCentroid nearest_centroid
        = { particle_metadata.current_cluster.idx, new_distance_2_to_current_cluster };

    // For each centroid, consider if it is closer than the current centroid.
    for (ui32 cluster_subset_idx = 0;
         cluster_subset_idx < Options.centroid_subset.k_prime;
         ++cluster_subset_idx)
    {
        ui32 cluster_idx = cluster_subset.indices[cluster_subset_idx];

        if (cluster_idx == particle_metadata.current_cluster.idx) continue;

        NBS_PRECISION centroid_distance_2
            = math::distance2(particle - clusters[cluster_idx].centroid);

        if (centroid_distance_2 < nearest_centroid.distance) {
            nearest_centroid.idx      = cluster_idx;
            nearest_centroid.distance = centroid_distance_2;
        }
    }

    return nearest_centroid;
}

template <nbs::Particle ParticleType, nbs::cluster::KMeansOptions Options>
nbs::cluster::detail::NearestCentroid
nbs::cluster::detail::nearest_centroid_and_build_list(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters,
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
        return math::distance2(particle - clusters[idx].centroid);
    };

    // Sort indices according to distance to particle.
    std::ranges::sort(
        buffers.nearest_centroid_indices, std::less<>{}, index_to_distance
    );

    // Set cluster subset.
    for (ui32 idx = 0; idx < Options.centroid_subset.k_prime; ++idx) {
        cluster_subset.indices[idx] = buffers.nearest_centroid_indices[idx];
    }

    return NearestCentroid{ buffers.nearest_centroid_indices[0],
                            math::distance2(
                                particle
                                - clusters[buffers.nearest_centroid_indices[0]].centroid
                            ) };
}

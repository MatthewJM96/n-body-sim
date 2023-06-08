template <nbs::Particle ParticleType, nbs::cluster::KMeansOptions Options>
nbs::cluster::detail::NearestCentroid nbs::cluster::detail::nearest_centroid(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters
) {
    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    NBS_PRECISION new_distance_2_to_current_cluster = math::distance2(
        particle - clusters[particle_metadata.current_cluster.idx].centroid
    );
    if (new_distance_2_to_current_cluster < particle_metadata.current_cluster.distance)
    {
        return NearestCentroid{ particle_metadata.current_cluster.idx,
                                new_distance_2_to_current_cluster };
    }

    NearestCentroid nearest_centroid
        = { particle_metadata.current_cluster.idx, new_distance_2_to_current_cluster };

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
    // Optimisation by early back out of search if previous nearest centroid has got
    // closer - in which case it is guaranteed to still be the nearest centroid.
    //     This is based on the paper "An Efficient Enhanced k-means Clustering
    //     Algorithm" by Fahim A.M., Salem A.M., Torkey F.A., and Ramadan M.A.
    NBS_PRECISION new_distance_2_to_current_cluster = math::distance2(
        particle - clusters[particle_metadata.current_cluster.idx].centroid
    );
    if (new_distance_2_to_current_cluster < particle_metadata.current_cluster.distance)
    {
        return NearestCentroid{ particle_metadata.current_cluster.idx,
                                new_distance_2_to_current_cluster };
    }

    NearestCentroid nearest_centroid
        = { particle_metadata.current_cluster.idx, new_distance_2_to_current_cluster };

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

#include <map>

template <nbs::Particle ParticleType, nbs::cluster::KMeansOptions Options>
nbs::cluster::detail::NearestCentroid
nbs::cluster::detail::nearest_centroid_and_build_list(
    const ParticleType&            particle,
    const ParticleClusterMetadata& particle_metadata,
    const Cluster<ParticleType>*   clusters,
    OUT detail::NearestCentroidList& cluster_subset
) {
    // Optimisation by making a subset of centroids to consider for a given particle.
    //     This is based on the paper "Faster k-means Cluster Estimation" by Khandelwal
    //     S., Awekar A.
    std::multimap<NBS_PRECISION, ui32> centroid_list;

    NearestCentroid nearest_centroid{
        particle_metadata.current_cluster.idx,
        math::distance2(
            particle - clusters[particle_metadata.current_cluster.idx].centroid
        )
    };

    for (ui32 cluster_idx = 0; cluster_idx < Options.cluster_count; ++cluster_idx) {
        if (cluster_idx == particle_metadata.current_cluster.idx) continue;

        NBS_PRECISION centroid_distance_2
            = math::distance2(particle - clusters[cluster_idx].centroid);

        centroid_list.insert(std::make_pair(centroid_distance_2, cluster_idx));

        if (centroid_distance_2 < nearest_centroid.distance) {
            nearest_centroid.idx      = cluster_idx;
            nearest_centroid.distance = centroid_distance_2;
        }
    }

    auto it = centroid_list.begin();
    for (ui32 idx = 0; idx < Options.centroid_subset.k_prime; ++idx) {
        cluster_subset.indices[idx] = *it++;
    }

    return nearest_centroid;
}

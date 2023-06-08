template <nbs::cluster::KMeansOptions Options>
void nbs::cluster::allocate_kmeans_buffers(
    OUT CALLER_DELETE KMeansBuffers<Options>& buffers
) {
    buffers.particle_metadata
        = new detail::ParticleClusterMetadata[Options.particle_count];
    buffers.cluster_modified_in_iteration = new bool[Options.cluster_count];

    if constexpr (Options.centroid_subset_optimisation) {
        buffers.nearest_centroid_indices = new size_t[Options.cluster_count];
        buffers.nearest_centroid_lists
            = new detail::NearestCentroidList[Options.particle_count];

        for (size_t i = 0; i < Options.particle_count; ++i) {
            buffers.nearest_centroid_lists[i]
                = new ui32[Options.centroid_subset.k_prime];
        }
    }
}

template <nbs::cluster::KMeansOptions Options>
void nbs::cluster::deallocate_kmeans_buffers(
    OUT CALLER_DELETE KMeansBuffers<Options> buffers
) {
    if constexpr (Options.centroid_subset_optimisation) {
        for (size_t i = 0; i < Options.particle_count; ++i) {
            delete[] buffers.nearest_centroid_lists[i];
        }

        delete[] buffers.nearest_centroid_lists;
        delete[] buffers.nearest_centroid_indices;
    }

    delete[] buffers.cluster_modified_in_iteration;
    delete[] buffers.particle_metadata;
}

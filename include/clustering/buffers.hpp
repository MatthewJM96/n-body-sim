#ifndef N_BODY_SIM_CLUSTERING_BUFFERS_HPP
#define N_BODY_SIM_CLUSTERING_BUFFERS_HPP

#pragma once

#include "clustering/options.hpp"

namespace nbs {
    namespace cluster {
        namespace detail {
            struct NearestCentroid {
                ui32          idx;
                NBS_PRECISION distance;
            };

            struct NearestCentroidList {
                ui32* indices;
            };
        }  // namespace detail

        template <KMeansOptions, typename = void>
        struct KMeansBuffers;

        template <KMeansOptions Options>
        struct KMeansBuffers<
            Options,
            typename std::enable_if_t<!Options.centroid_subset_optimisation>> {
            detail::NearestCentroid* particle_nearest_centroid;
            bool*                    cluster_modified_in_iteration;
        };

        template <KMeansOptions Options>
        struct KMeansBuffers<
            Options,
            typename std::enable_if_t<Options.centroid_subset_optimisation>> {
            detail::NearestCentroid*     particle_nearest_centroid;
            bool*                        cluster_modified_in_iteration;
            detail::NearestCentroidList* nearest_centroid_lists;
            ui32*                        nearest_centroid_indices;
        };

        template <KMeansOptions Options>
        void allocate_kmeans_buffers(OUT CALLER_DELETE KMeansBuffers<Options> buffers);

        template <KMeansOptions Options>
        void deallocate_kmeans_buffers(OUT CALLER_DELETE KMeansBuffers<Options> buffers
        );
    }  // namespace cluster
}  // namespace nbs

#include "buffers.inl"

#endif  // N_BODY_SIM_CLUSTERING_BUFFERS_HPP

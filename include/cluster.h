#ifndef N_BODY_SIM_CLUSTER_H
#define N_BODY_SIM_CLUSTER_H

#pragma once

namespace nbs {
    namespace cluster {
        template <
            typename Precision,
            typename =
                typename std::enable_if<std::is_floating_point<Precision>::value>::type>
        struct Member {
            Precision x, y, z;
        };

        template <
            typename Precision,
            typename =
                typename std::enable_if<std::is_floating_point<Precision>::value>::type>
        struct Cluster {
            Member<Precision>  centroid;
            Member<Precision>* members;
            ui32               member_count;
        };

        struct KMeansOptions {
            ui32 max_iterations                       = 100;
            ui32 acceptable_changes_per_iteration     = 0;
            bool front_loaded                         = false;
            bool no_approaching_centroid_optimisation = false;
            bool no_centroid_subset_optimisation      = true;

            struct {
                ui32 k_prime = 30;
            } centroid_subset;
        };

        template <typename Precision>
        Precision
        member_distance_2(const Member<Precision>& lhs, const Member<Precision>& rhs);

        template <typename Precision>
        void
        kpp(Member<Precision>* members,
            ui32               member_count,
            Member<Precision>* centroids,
            ui32               centroid_count);

        template <typename Precision>
        void k_means(
            const Cluster<Precision>* initial_clusters,
            ui32                      initial_cluster_count,
            const KMeansOptions&      options,
            OUT Cluster<Precision>*& clusters
        );

        namespace impl {
            template <
                typename Precision,
                typename = typename std::enable_if<
                    std::is_floating_point<Precision>::value>::type>
            struct NearestCentroid {
                ui32      idx;
                Precision distance;
            };

            template <
                typename Precision,
                typename = typename std::enable_if<
                    std::is_floating_point<Precision>::value>::type>
            struct MemberClusterMetadata {
                ui32 initial_cluster_idx;
                ui32 initial_member_idx;

                NearestCentroid<Precision> current_cluster;
            };

            struct NearestCentroidList {
                ui32* indices;
            };

            template <typename Precision>
            struct NearestCentroidAndList {
                NearestCentroid<Precision> centroid;
                NearestCentroidList        list;
            };

            template <typename Precision>
            NearestCentroid<Precision> nearest_centroid(
                const Member<Precision>&                member,
                const MemberClusterMetadata<Precision>& member_metadata,
                const Cluster<Precision>*               clusters,
                ui32                                    cluster_count
            );
            template <typename Precision>
            NearestCentroid<Precision> nearest_centroid_from_subset(
                const Member<Precision>&                member,
                const MemberClusterMetadata<Precision>& member_metadata,
                const Cluster<Precision>*               clusters,
                NearestCentroidList                     cluster_subset,
                ui32                                    cluster_count
            );

            template <typename Precision>
            NearestCentroidAndList<Precision> nearest_centroid_and_build_list(
                const Member<Precision>&                member,
                const MemberClusterMetadata<Precision>& member_metadata,
                const Cluster<Precision>*               clusters,
                ui32                                    cluster_count,
                ui32                                    subset_count
            );
        }  // namespace impl
    }      // namespace cluster
}  // namespace nbs

#include "cluster.inl"

#endif  // N_BODY_SIM_CLUSTER_H

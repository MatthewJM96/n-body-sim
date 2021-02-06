#ifndef N_BODY_SIM_CLUSTER_H
#define N_BODY_SIM_CLUSTER_H

#pragma once

#include <type_traits>

#include "types.h"

namespace cluster {
    template <typename Precision, typename = typename std::enable_if<std::is_floating_point<Precision>::value>::type>
    struct Member {
        Precision x, y, z;
    };

    template <typename Precision, typename = typename std::enable_if<std::is_floating_point<Precision>::value>::type>
    struct Cluster {
        Member<Precision> centroid;
        Member<Precision>* members;
        ui32 member_count;
    };

    struct KMeansOptions {
        ui32 max_iterations;
        ui32 acceptable_changes_per_iteration;
    };

    template <typename Precision>
    Precision member_distance_2(const Member<Precision>& lhs, const Member<Precision>& rhs);

    template <typename Precision>
    void kpp(Member<Precision>* members, ui32 member_count, Member<Precision>* centroids, ui32 centroid_count);

    template <typename Precision>
    void k_means(const Cluster<Precision>* initial_clusters, ui32 initial_cluster_count, const KMeansOptions& options, OUT Cluster<Precision>*& clusters, bool front_loaded = false);

    namespace impl {
        template <typename Precision>
        ui32 nearest_centroid(const Member<Precision>& member, const Cluster<Precision>* clusters, ui32 cluster_count);

        struct MemberClusterMetadata {
            ui32 initial_cluster_idx;
            ui32 initial_member_idx;
            ui32 current_cluster_idx;
        };
    };
};

#include "cluster.inl"

#endif // N_BODY_SIM_CLUSTER_H

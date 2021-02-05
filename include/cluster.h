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
        Member<Precision>* members;
        ui32 member_count;
    }

    template <typename Precision>
    Precision member_distance_2(Member<Precision> lhs, Member<Precision> rhs);

    template <typename Precision>
    void kpp(Member<Precision>* members, ui32 member_count, Member<Precision>* centroids, ui32 centroid_count);

    template <typename Precision>
    void k_means(Member<Precision>* members, ui32 member_count, Member<Precision>* centroids, Cluster<Precision>* clusters, ui32 centroid_count);
};

#include "cluster.inl"

#endif // N_BODY_SIM_CLUSTER_H

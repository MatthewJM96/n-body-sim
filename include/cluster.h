#ifndef N_BODY_SIM_CLUSTER_H
#define N_BODY_SIM_CLUSTER_H

#pragma once

#include <type_traits>

#include "types.h"

namespace cluster {
    template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
    struct Member {
        T x, y, z;
    };

    template <typename Precision>
    Precision member_distance_2(Member<Precision> lhs, Member<Precision> rhs);

    template <typename Precision>
    void cluster::kpp(Member<Precision>* members, ui32 member_count, Member<Precision>* centroids, ui32 centroid_count);
};

#include "cluster.inl"

#endif // N_BODY_SIM_CLUSTER_H

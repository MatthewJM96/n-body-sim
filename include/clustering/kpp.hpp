#ifndef N_BODY_SIM_CLUSTERING_KPP_HPP
#define N_BODY_SIM_CLUSTERING_KPP_HPP

#pragma once

#include "particle.hpp"

namespace nbs {
    namespace cluster {
        template <ClusteredParticle ParticleType, KMeansOptions Options>
        void kpp(const ParticleType* particles, IN OUT Cluster<ParticleType>* clusters);
    }  // namespace cluster
}  // namespace nbs

#include "kpp.inl"

#endif  // N_BODY_SIM_CLUSTERING_KPP_HPP

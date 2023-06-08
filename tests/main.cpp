#include "stdafx.h"

#include "clustering/clustering.hpp"

using namespace nbs;

struct MyParticle {
    f32v3  position;
    size_t cluster_metadata_idx;
};

template <size_t ParticleCount, size_t ClusterCount>
void do_a_cluster_job() {
    constexpr cluster::KMeansOptions options
        = { .particle_count = ParticleCount, .cluster_count = ClusterCount };

    // Allocate particles.
    MyParticle* particles = new MyParticle[ParticleCount];

    // Set up particles.
    std::default_random_engine          generator;
    std::uniform_real_distribution<f32> distribution(-1000.0f, 1000.0f);
    for (size_t i = 0; i < ParticleCount; ++i) {
        particles[i].cluster_metadata_idx = i;
        particles[i].position             = f32v3(
            distribution(generator), distribution(generator), distribution(generator)
        );
    }

    // Allocate clusters.
    cluster::Cluster<MyParticle>* clusters
        = new cluster::Cluster<MyParticle>[ClusterCount * 2];

    // Do kpp initialisation.
    cluster::kpp<MyParticle, options>(particles, clusters);

    // Quick check.
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << clusters[i].centroid.position.x << " - "
                  << clusters[i].centroid.position.y << " - "
                  << clusters[i].centroid.position.z << std::endl;
    }
}

int main() {
    do_a_cluster_job<1000, 10>();
}

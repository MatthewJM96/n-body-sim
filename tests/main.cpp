#include "stdafx.h"

#include "clustering/clustering.hpp"

using namespace nbs;

struct MyParticle2D {
    f32v2  position;
    size_t cluster_metadata_idx;
};

struct MyParticle {
    f32v3  position;
    size_t cluster_metadata_idx;
};

template <size_t ParticleCount, size_t ClusterCount>
void do_a_cluster_job_dim_2() {
    constexpr cluster::KMeansOptions options = { .particle_count = ParticleCount,
                                                 .cluster_count  = ClusterCount,
                                                 .front_loaded   = true };

    // Allocate particles.
    MyParticle2D* particles = new MyParticle2D[ParticleCount];

    // Set up particles.
    std::default_random_engine          generator;
    std::uniform_real_distribution<f32> distribution(-1000.0f, 1000.0f);
    for (size_t i = 0; i < ParticleCount; ++i) {
        particles[i].cluster_metadata_idx = i;
        particles[i].position = f32v2(distribution(generator), distribution(generator));
    }

    // Allocate clusters.
    cluster::Cluster<2, MyParticle2D>* clusters
        = new cluster::Cluster<2, MyParticle2D>[ClusterCount * 2];

    // Do kpp initialisation.
    cluster::kpp<2, MyParticle2D, options>(particles, clusters);

    // Quick check.
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << clusters[i].centroid.position.x << " - "
                  << clusters[i].centroid.position.y << std::endl;
    }

    // Front load into first cluster.
    clusters[0].particle_count  = ParticleCount;
    clusters[0].particle_offset = 0;

    // Allocate buffers used for K-means.
    cluster::KMeansBuffers<options> buffers;
    cluster::allocate_kmeans_buffers<options>(buffers);

    // Do k_means.
    cluster::k_means<2, MyParticle2D, options>(
        particles, clusters, clusters + ClusterCount, buffers
    );

    // Quick check.
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << clusters[i + ClusterCount].centroid.position.x << " - "
                  << clusters[i + ClusterCount].centroid.position.y << std::endl;
    }
}

template <size_t ParticleCount, size_t ClusterCount>
void do_a_cluster_job_dim_3() {
    constexpr cluster::KMeansOptions options = { .particle_count = ParticleCount,
                                                 .cluster_count  = ClusterCount,
                                                 .front_loaded   = true };

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
    cluster::Cluster<3, MyParticle>* clusters
        = new cluster::Cluster<3, MyParticle>[ClusterCount * 2];

    // Do kpp initialisation.
    cluster::kpp<3, MyParticle, options>(particles, clusters);

    // Quick check.
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << clusters[i].centroid.position.x << " - "
                  << clusters[i].centroid.position.y << " - "
                  << clusters[i].centroid.position.z << std::endl;
    }

    // Front load into first cluster.
    clusters[0].particle_count  = ParticleCount;
    clusters[0].particle_offset = 0;

    // Allocate buffers used for K-means.
    cluster::KMeansBuffers<options> buffers;
    cluster::allocate_kmeans_buffers<options>(buffers);

    // Do k_means.
    cluster::k_means<3, MyParticle, options>(
        particles, clusters, clusters + ClusterCount, buffers
    );

    // Quick check.
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << clusters[i + ClusterCount].centroid.position.x << " - "
                  << clusters[i + ClusterCount].centroid.position.y << " - "
                  << clusters[i + ClusterCount].centroid.position.z << std::endl;
    }
}

int main() {
    do_a_cluster_job_dim_2<1000, 10>();

    do_a_cluster_job_dim_3<1000, 10>();
}

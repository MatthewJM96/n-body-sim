#include "stdafx.h"

#include <chrono>

#include "clustering/clustering.hpp"

#include "2D_clustering_viewer.hpp"
#include "my_particles.hpp"

#include "data/a1.hpp"

#include "statistics/average_cluster_distance.hpp"

using namespace nbs;

// TODO(Matthew): Make timing more robust.
// TODO(Matthew): Implement a test case of clustering, then moving particles, and then
//                timing the reclustering after this. Parameterise the degree of
//                particle movement and get a vibe of the consequence of more or less
//                movement.
// TODO(Matthew): Implement larger test cases and trial optimisations.
// TODO(Matthew): Implement retrialing kpp over X attempts to find best seeding.
//                  Do the same with kpp + k_means combined? Ultimately for simulation
//                  our first centroids must be GOOD.
// TODO(Matthew): Validate that particle motion doesn't sufficiently screw up the best
//                centroids that the time to get them isn't worth it.
// TODO(Matthew): Write a first-pass force update for particles.
// TODO(Matthew): Write one or two nice distributions for particles.

template <size_t ClusterCount, size_t Iterations>
void do_a_cluster_job_a1(
    MyParticle2D*& particles, cluster::Cluster<2, MyParticle2D>*& clusters
) {
    constexpr cluster::KMeansOptions options
        = { .particle_count                    = 7500,
            .cluster_count                     = ClusterCount,
            .max_iterations                    = 100,
            .front_loaded                      = true,
            .approaching_centroid_optimisation = false };

    // Allocate particles.
    particles = new MyParticle2D[7500];

    // Set up particles.
    for (size_t i = 0; i < 7500; ++i) {
        particles[i].cluster_metadata_idx = i;
        particles[i].position             = A1_DATA[i];
    }

    // Allocate clusters.
    clusters = new cluster::Cluster<2, MyParticle2D>[ClusterCount * 2];

    nbs::i64 total_us = 0;
    for (int iteration = 0; iteration < Iterations; ++iteration) {
        // Do kpp initialisation.
        cluster::kpp<2, MyParticle2D, options>(particles, clusters);

        // // Quick check.
        // std::cout << "    kpp centroids:" << std::endl;
        // // clang-format off
        // for (size_t i = 0; i < ClusterCount; ++i) {
        //     std::cout << "        " << clusters[i].centroid.position.x << " - "
        //                             << clusters[i].centroid.position.y << std::endl;
        // }
        // // clang-format on

        // Front load into first cluster.
        clusters[0].particle_count  = 7500;
        clusters[0].particle_offset = 0;

        // Allocate buffers used for K-means.
        cluster::KMeansBuffers<options> buffers;
        cluster::allocate_kmeans_buffers<options>(buffers);

        auto start = std::chrono::high_resolution_clock::now();
        // Do k_means.
        cluster::k_means<2, MyParticle2D, options>(
            particles, clusters, clusters + ClusterCount, buffers
        );
        auto duration = std::chrono::high_resolution_clock::now() - start;
        total_us
            += std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    }

    if constexpr (Iterations != 1) {
        std::cout << "Average time to cluster: "
                  << static_cast<f32>(total_us) / static_cast<f32>(Iterations) << "us"
                  << std::endl;
    }

    std::cout
        << "Average particle distance to cluster: "
        << statistics::
               calculate_average_cluster_distance<2, MyParticle2D, 7500, ClusterCount>(
                   particles, clusters
               )
        << std::endl;

    // // Quick check.
    // std::cout << "    k_means centroids:" << std::endl;
    // // clang-format off
    // for (size_t i = 0; i < ClusterCount; ++i) {
    //     std::cout << "        " << clusters[i + ClusterCount].centroid.position.x <<
    //     " - "
    //                             << clusters[i + ClusterCount].centroid.position.y <<
    //                             std::endl;
    // }
    // // clang-format on
}

template <size_t ParticleCount, size_t ClusterCount>
void do_a_cluster_job_dim_2(
    MyParticle2D*& particles, cluster::Cluster<2, MyParticle2D>*& clusters
) {
    constexpr cluster::KMeansOptions options = { .particle_count = ParticleCount,
                                                 .cluster_count  = ClusterCount,
                                                 .max_iterations = 100,
                                                 .front_loaded   = true };

    // Allocate particles.
    particles = new MyParticle2D[ParticleCount];

    // Set up particles.
    std::default_random_engine          generator;
    std::uniform_real_distribution<f32> distribution(-1000.0f, 1000.0f);
    for (size_t i = 0; i < ParticleCount; ++i) {
        particles[i].cluster_metadata_idx = i;
        particles[i].position = f32v2(distribution(generator), distribution(generator));
    }

    // Allocate clusters.
    clusters = new cluster::Cluster<2, MyParticle2D>[ClusterCount * 2];

    // Do kpp initialisation.
    cluster::kpp<2, MyParticle2D, options>(particles, clusters);

    // Quick check.
    std::cout << "    kpp centroids:" << std::endl;
    // clang-format off
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << "        " << clusters[i].centroid.position.x << " - "
                                << clusters[i].centroid.position.y << std::endl;
    }
    // clang-format on

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
    std::cout << "    k_means centroids:" << std::endl;
    // clang-format off
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << "        " << clusters[i + ClusterCount].centroid.position.x << " - "
                                << clusters[i + ClusterCount].centroid.position.y << std::endl;
    }
    // clang-format on
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
    std::cout << "    kpp centroids:" << std::endl;
    // clang-format off
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << "        " << clusters[i].centroid.position.x << " - "
                                << clusters[i].centroid.position.y << " - "
                                << clusters[i].centroid.position.z << std::endl;
    }
    // clang-format on

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
    // clang-format off
    std::cout << "    k_means centroids:" << std::endl;
    for (size_t i = 0; i < ClusterCount; ++i) {
        std::cout << "        " << clusters[i + ClusterCount].centroid.position.x << " - "
                                << clusters[i + ClusterCount].centroid.position.y << " - "
                                << clusters[i + ClusterCount].centroid.position.z << std::endl;
    }
    // clang-format on
}

void do_2D_uniform_distribution_case() {
#define PARTICLE_COUNT 1000
#define CLUSTER_COUNT  10

    MyParticle2D*                      particles;
    cluster::Cluster<2, MyParticle2D>* clusters;

    do_a_cluster_job_dim_2<PARTICLE_COUNT, CLUSTER_COUNT>(particles, clusters);

#undef PARTICLE_COUNT
#undef CLUSTER_COUNT
}

void do_3d_uniform_distribution_case() {
#define PARTICLE_COUNT 1000
#define CLUSTER_COUNT  10

    do_a_cluster_job_dim_3<PARTICLE_COUNT, CLUSTER_COUNT>();

#undef PARTICLE_COUNT
#undef CLUSTER_COUNT
}

void do_a1_dataset_case() {
    MyParticle2D*                      particles;
    cluster::Cluster<2, MyParticle2D>* clusters;

    do_a_cluster_job_a1<50, 1>(particles, clusters);

    make_2d_cluster_view(particles, clusters);
}

void do_a1_dataset_performance_case() {
    MyParticle2D*                      particles;
    cluster::Cluster<2, MyParticle2D>* clusters;

    do_a_cluster_job_a1<50, 1000>(particles, clusters);

    make_2d_cluster_view(particles, clusters);
}

int main() {
    std::cout << "N-Body Simulator Menu:\n"
                 "  - 2D Uniform Distribution Case (1)\n"
                 "  - 3D Uniform Distribution Case (2)\n"
                 "  - A1 Dataset Case              (3)\n"
                 "  - A1 Dataset Performance Case  (4)\n"
              << std::endl;

    char resp;
    std::cin >> resp;

    if (resp == '1') {
        do_2D_uniform_distribution_case();
    } else if (resp == '2') {
        do_3d_uniform_distribution_case();
    } else if (resp == '3') {
        do_a1_dataset_case();
    } else if (resp == '4') {
        do_a1_dataset_performance_case();
    }
}

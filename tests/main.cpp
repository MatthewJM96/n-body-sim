#include "stdafx.h"

#include <chrono>

#include "clustering/clustering.hpp"

#include "2D_clustering_viewer.hpp"
#include "my_particles.hpp"

#include "data/a1.hpp"

#include "statistics/average_cluster_distance.hpp"

#include "forces/gravity.hpp"

using namespace nbs;

// TODO(Matthew): Make timing more robust.
// TODO(Matthew): Implement a test case of clustering, then moving particles, and then
//                timing the reclustering after this. Parameterise the degree of
//                particle movement and get a vibe of the consequence of more or less
//                movement.
// TODO(Matthew): Implement larger test cases and trial optimisations.
// TODO(Matthew): Validate that particle motion doesn't sufficiently screw up the best
//                centroids that the time to get them isn't worth it.
// TODO(Matthew): Write a first-pass force update for particles.
// TODO(Matthew): Write one or two nice distributions for particles.
// TODO(Matthew): Explore bisecting K-means (K-means for 2 centroids replacing each
//                centroid in a current iteration, starting with 1 centroid over whole
//                dataset).

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
                   particles, clusters + ClusterCount
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

template <size_t ClusterCount, size_t Attempts>
void do_optimise_kpp_a1_job(
    MyParticle2D*& particles, cluster::Cluster<2, MyParticle2D>*& clusters
) {
    constexpr cluster::KMeansOptions options
        = { .particle_count                    = 7500,
            .cluster_count                     = ClusterCount,
            .max_iterations                    = 100,
            .front_loaded                      = true,
            .approaching_centroid_optimisation = false };

    std::random_device rand_dev;

    // Allocate particles.
    particles = new MyParticle2D[7500];

    // Allocate clusters.
    clusters = new cluster::Cluster<2, MyParticle2D>[ClusterCount * 2];

    f32  current_best_avg_dist = std::numeric_limits<f32>::max();
    ui32 current_best_seed     = 0;

    nbs::i64 total_us = 0;
    for (int iteration = 0; iteration < Attempts; ++iteration) {
        ui32 seed = rand_dev();

        // Set up particles.
        for (size_t i = 0; i < 7500; ++i) {
            particles[i].cluster_metadata_idx = i;
            particles[i].position             = A1_DATA[i];
        }

        // Do kpp initialisation.
        cluster::kpp<2, MyParticle2D, options>(particles, clusters, &seed);

        // Front load into first cluster.
        clusters[0].particle_count  = 7500;
        clusters[0].particle_offset = 0;

        // Allocate buffers used for K-means.
        cluster::KMeansBuffers<options> buffers;
        cluster::allocate_kmeans_buffers<options>(buffers);

        // Do k_means.
        cluster::k_means<2, MyParticle2D, options>(
            particles, clusters, clusters + ClusterCount, buffers
        );

        f32 avg_dist = statistics::
            calculate_average_cluster_distance<2, MyParticle2D, 7500, ClusterCount>(
                particles, clusters + ClusterCount
            );
        if (avg_dist < current_best_avg_dist) {
            current_best_avg_dist = avg_dist;
            current_best_seed     = seed;
        }
    }

    // Now do best again for visualisation.

    // Set up particles.
    for (size_t i = 0; i < 7500; ++i) {
        particles[i].cluster_metadata_idx = i;
        particles[i].position             = A1_DATA[i];
    }

    // Do kpp initialisation.
    cluster::kpp<2, MyParticle2D, options>(particles, clusters, &current_best_seed);

    // Front load into first cluster.
    clusters[0].particle_count  = 7500;
    clusters[0].particle_offset = 0;

    // Allocate buffers used for K-means.
    cluster::KMeansBuffers<options> buffers;
    cluster::allocate_kmeans_buffers<options>(buffers);

    // Do k_means.
    cluster::k_means<2, MyParticle2D, options>(
        particles, clusters, clusters + ClusterCount, buffers
    );

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

    std::cout << "Achieved average particle distance to cluster: "
              << current_best_avg_dist << std::endl;
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

template <size_t ClusterCount>
void do_run_sim_step(
    MyParticle2D* particles, cluster::Cluster<2, MyParticle2D>* clusters
) {
    for (size_t cluster_idx = 0; cluster_idx < ClusterCount; ++cluster_idx) {
        const auto& cluster = clusters[cluster_idx];

        for (size_t offset = 0; offset < cluster.particle_count; ++offset) {
            auto& particle = particles[cluster.particle_offset + offset];

            particle.force = {};
        }

        for (size_t p1_offset = 0; p1_offset < cluster.particle_count; ++p1_offset) {
            auto& particle_1 = particles[cluster.particle_offset + p1_offset];
            for (size_t p2_offset = p1_offset + 1; p2_offset < cluster.particle_count;
                 ++p2_offset)
            {
                auto& particle_2 = particles[cluster.particle_offset + p2_offset];

                f32 distance_2
                    = math::distance2(particle_1.position, particle_2.position);

                // f32 force = forces::grav_with_repulsion_6<1000>(distance_2);
                f32 force = forces::grav(distance_2);

                particle_1.force
                    += math::normalize(particle_2.position - particle_1.position)
                       * force;
                particle_2.force
                    += math::normalize(particle_1.position - particle_2.position)
                       * force;
            }

            for (size_t other_cluster_idx = 0; other_cluster_idx < ClusterCount;
                 ++other_cluster_idx)
            {
                const auto& other_cluster = clusters[other_cluster_idx];

                f32 distance_2
                    = math::distance2(particle_1.position, cluster.centroid.position);

                // f32 force = forces::grav_with_repulsion_6<1000>(distance_2);
                f32 force = forces::grav(distance_2);

                particle_1.force
                    += math::normalize(cluster.centroid.position - particle_1.position)
                       * force * static_cast<f32>(other_cluster.particle_count);
            }
        }

        for (size_t offset = 0; offset < cluster.particle_count; ++offset) {
            auto& particle = particles[cluster.particle_offset + offset];

            const f32 t_fact = 100.0f;

            particle.velocity += particle.force * t_fact;
            particle.position
                += particle.velocity * t_fact - 0.5f * particle.force * t_fact * t_fact;
        }
    }
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

void do_a1_dataset_optimise_kpp_case() {
    const f32v4 clip_rect = f32v4(-1000.0f, 65000.0f, -1000.0f, 66000.0f);

    MyParticle2D*                      particles;
    cluster::Cluster<2, MyParticle2D>* clusters;

    do_optimise_kpp_a1_job<50, 100>(particles, clusters);

    make_2d_cluster_view(particles, clusters, &clip_rect);

    for (size_t i = 0; i < 7500; ++i) {
        particles[i].velocity = {};
        particles[i].force    = {};
    }

    for (size_t i = 0; i < 10; ++i) {
        do_run_sim_step<50>(particles, clusters + 50);
    }

    make_2d_cluster_view(particles, clusters, &clip_rect);

    for (size_t i = 0; i < 20; ++i) {
        do_run_sim_step<50>(particles, clusters + 50);
    }

    make_2d_cluster_view(particles, clusters, &clip_rect);

    for (size_t i = 0; i < 50; ++i) {
        do_run_sim_step<50>(particles, clusters + 50);
    }

    make_2d_cluster_view(particles, clusters, &clip_rect);

    for (size_t i = 0; i < 80; ++i) {
        do_run_sim_step<50>(particles, clusters + 50);
    }

    make_2d_cluster_view(particles, clusters, &clip_rect);
}

int main() {
    std::cout << "N-Body Simulator Menu:\n"
                 "  - 2D Uniform Distribution Case (1)\n"
                 "  - 3D Uniform Distribution Case (2)\n"
                 "  - A1 Dataset Case              (3)\n"
                 "  - A1 Dataset Performance Case  (4)\n"
                 "  - A1 Dataset Optimise KPP Case (5)\n"
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
    } else if (resp == '5') {
        do_a1_dataset_optimise_kpp_case();
    }
}

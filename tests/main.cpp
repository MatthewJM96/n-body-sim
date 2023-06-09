#include "stdafx.h"

#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <implot.h>

#include "clustering/clustering.hpp"

#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

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
void do_a_cluster_job_dim_2(
    MyParticle2D*& particles, cluster::Cluster<2, MyParticle2D>*& clusters
) {
    constexpr cluster::KMeansOptions options = { .particle_count = ParticleCount,
                                                 .cluster_count  = ClusterCount,
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

int main() {
    MyParticle2D*                      particles;
    cluster::Cluster<2, MyParticle2D>* clusters;

#define PARTICLE_COUNT 1000
#define CLUSTER_COUNT  10

    // std::cout << glGetString(GL_VERSION) << std::endl;

    std::cout << "2D case:" << std::endl;
    do_a_cluster_job_dim_2<PARTICLE_COUNT, CLUSTER_COUNT>(particles, clusters);

    std::cout << std::endl << std::endl;

    std::cout << "3D case:" << std::endl;
    do_a_cluster_job_dim_3<PARTICLE_COUNT, CLUSTER_COUNT>();

    glfwInit();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(1000, 800, "My Window", nullptr, nullptr);

    glfwMakeContextCurrent(window);

    IMGUI_CHECKVERSION();

    auto gui_ctx = ImGui::CreateContext();

    ImGui_ImplGlfw_InitForOpenGL(window, false);
    ImGui_ImplOpenGL3_Init();

    auto plt_ctx = ImPlot::CreateContext();

    ImGui::SetCurrentContext(gui_ctx);
    ImPlot::SetCurrentContext(plt_ctx);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowSize(ImVec2(1600, 1200));
        if (ImGui::Begin("My SubWindow")) {
            if (ImPlot::BeginPlot("My Plot", ImVec2(1500, 1100))) {
                ImPlot::PlotScatterG(
                    "particles",
                    [](int idx, void* data) {
                        auto particle = reinterpret_cast<MyParticle2D*>(data)[idx];
                        ImPlotPoint p = { particle.position.x, particle.position.y };
                        return p;
                    },
                    reinterpret_cast<void*>(particles),
                    PARTICLE_COUNT
                );

                ImPlot::SetNextMarkerStyle(ImPlotMarker_Diamond, 6.0f);
                ImPlot::PlotScatterG(
                    "centroids",
                    [](int idx, void* data) {
                        auto cluster
                            = reinterpret_cast<cluster::Cluster<2, MyParticle2D>*>(data
                            )[idx];
                        ImPlotPoint p = { cluster.centroid.position.x,
                                          cluster.centroid.position.y };
                        return p;
                    },
                    reinterpret_cast<void*>(clusters + CLUSTER_COUNT),
                    CLUSTER_COUNT
                );

                ImPlot::EndPlot();
            }
            ImGui::End();
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glfwSwapBuffers(window);
    }

    ImPlot::DestroyContext(plt_ctx);
    ImGui::DestroyContext(gui_ctx);
}

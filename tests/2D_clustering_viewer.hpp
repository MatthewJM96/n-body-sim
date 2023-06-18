#ifndef N_BODY_SIM_TESTS_2D_CLUSTERING_VIEWER_HPP
#define N_BODY_SIM_TESTS_2D_CLUSTERING_VIEWER_HPP

#pragma once

#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <implot.h>

#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "clustering/clustering.hpp"

#include "my_particles.hpp"

void make_2d_cluster_view(
    MyParticle2D*                           particles,
    nbs::cluster::Cluster<2, MyParticle2D>* clusters,
    const nbs::f32v4*                       clip_rect = nullptr
) {
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
                if (clip_rect)
                    ImPlot::SetupAxesLimits(
                        clip_rect->x, clip_rect->y, clip_rect->z, clip_rect->w
                    );

                ImPlot::PlotScatterG(
                    "particles",
                    [](int idx, void* data) {
                        auto particle = reinterpret_cast<MyParticle2D*>(data)[idx];
                        ImPlotPoint p = { particle.position.x, particle.position.y };
                        return p;
                    },
                    reinterpret_cast<void*>(particles),
                    7500
                );

                ImPlot::SetNextMarkerStyle(ImPlotMarker_Diamond, 6.0f);
                ImPlot::PlotScatterG(
                    "centroids",
                    [](int idx, void* data) {
                        auto cluster
                            = reinterpret_cast<nbs::cluster::Cluster<2, MyParticle2D>*>(
                                data
                            )[idx];
                        ImPlotPoint p = { cluster.centroid.position.x,
                                          cluster.centroid.position.y };
                        return p;
                    },
                    reinterpret_cast<void*>(clusters + 50),
                    50
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

#endif  // N_BODY_SIM_TESTS_2D_CLUSTERING_VIEWER_HPP

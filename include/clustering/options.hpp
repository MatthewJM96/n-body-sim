#ifndef N_BODY_SIM_CLUSTERING_OPTIONS_HPP
#define N_BODY_SIM_CLUSTERING_OPTIONS_HPP

#pragma once

namespace nbs {
    namespace cluster {
        struct KMeansOptions {
            ui32 particle_count                    = 1000;
            ui32 cluster_count                     = 10;
            ui32 max_iterations                    = 100;
            ui32 acceptable_changes_per_iteration  = 0;
            bool front_loaded                      = false;
            bool approaching_centroid_optimisation = true;
            bool centroid_subset_optimisation      = false;

            struct {
                ui32 k_prime    = 30;
                bool do_rebuild = true;
            } centroid_subset = {};
        };
    }  // namespace cluster
}  // namespace nbs

#endif  // N_BODY_SIM_CLUSTERING_OPTIONS_HPP

#ifndef N_BODY_SIM_TESTS_MY_PARTICLES_HPP
#define N_BODY_SIM_TESTS_MY_PARTICLES_HPP

#pragma once

struct MyParticle2D {
    nbs::f32v2 position;
    size_t     cluster_metadata_idx;
    nbs::f32v2 force;
    nbs::f32v2 velocity;
};

struct MyParticle {
    nbs::f32v3 position;
    size_t     cluster_metadata_idx;
    nbs::f32v2 force;
};

#endif  // N_BODY_SIM_TESTS_MY_PARTICLES_HPP

#ifndef N_BODY_SIM_FORCES_GRAV_WITH_REPULSION_HPP
#define N_BODY_SIM_FORCES_GRAV_WITH_REPULSION_HPP

#pragma once

namespace nbs {
    namespace forces {
        inline NBS_PRECISION grav(NBS_PRECISION distance_2) {
            return -1.0 / distance_2;
        }

        template <size_t Tightness>
        inline NBS_PRECISION grav_with_repulsion_6(NBS_PRECISION distance_2) {
            const NBS_PRECISION tightness = static_cast<NBS_PRECISION>(Tightness);

            return (-tightness / distance_2 + 1 / math::pow(distance_2, 3))
                   / (math::pow(tightness, 1.5) * -0.384900179459);
        }
    }  // namespace forces
}  // namespace nbs

#endif  // N_BODY_SIM_FORCES_GRAV_WITH_REPULSION_HPP

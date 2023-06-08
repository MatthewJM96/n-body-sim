#include "particle.h"

template <
    typename Precision,
    typename = typename std::enable_if<std::is_floating_point<Precision>::value>::type>
using Volume = Particle<Precision>;

template <
    typename Precision,
    typename = typename std::enable_if<std::is_floating_point<Precision>::value>::type>
class IParticleGenerator {
public:
    /**
     * \brief Generates a distribution of particles with a target particle count.
     *
     * \param particles The buffer of particles to fill up to the target particle count.
     * \param target_particle_count The target particle count for the distribution
     * generation, this will not be exceeded. \param box The box in which to place the
     * distribution of particles.
     *
     * \return The number of particles actually placed into the buffer.
     */
    virtual ui32 generate(
        Particle<Precision>* particles,
        ui32                 target_particle_count,
        Volume<Precision>    box
    ) = 0;
};

template <
    typename Precision,
    typename = typename std::enable_if<std::is_floating_point<Precision>::value>::type>
class CubeParticleGenerator {
public:
    /**
     * \brief Generates a distribution of particles that puts them into a cube layout
     * with equal distance along each axis.
     *
     * \param particles The buffer of particles to fill up to the target particle count.
     * \param target_particle_count The target particle count for the distribution
     * generation, this will not be exceeded. \param box The box in which to place the
     * distribution of particles.
     *
     * \return The number of particles actually placed into the buffer.
     */
    virtual ui32 generate(
        Particle<Precision>* particles,
        ui32                 target_particle_count,
        Volume<Precision>    box
    ) override {
        ui32 target_particle_count_on_axis
            = std::floor(std::cbrtf((float)target_particle_count));

        for (ui32 x = 0; x < target_particle_count_on_axis; ++x) {
            for (ui32 y = 0; y < target_particle_count_on_axis; ++y) {
                for (ui32 z = 0; z < target_particle_count_on_axis; ++z) {
                    ui32 buffer_idx = x * target_particle_count_on_axis
                                          * target_particle_count_on_axis
                                      + y * target_particle_count_on_axis + z;

                    if (buffer_idx + 1 >= target_particle_count)
                        return target_particle_count;

                    float x_ratio = (float)x / (float)target_particle_count_on_axis;
                    float y_ratio = (float)y / (float)target_particle_count_on_axis;
                    float z_ratio = (float)z / (float)target_particle_count_on_axis;

                    particles[buffer_idx] = Particle<Precision>{ x_ratio * box.x,
                                                                 y_ratio * box.y,
                                                                 z_ratio * box.z };
                }
            }
        }

        return target_particle_count;
    }
};

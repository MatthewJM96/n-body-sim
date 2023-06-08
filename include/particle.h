template <
    typename Precision,
    typename = typename std::enable_if<std::is_floating_point<Precision>::value>::type>
struct Particle {
    Precision x, y, z;
};

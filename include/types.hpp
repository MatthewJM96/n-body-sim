#ifndef N_BODY_SIM_TYPES_HPP
#define N_BODY_SIM_TYPES_HPP

#pragma once

namespace nbs {
    // Integral types.
    using i8   = int8_t;    ///< 8-bit signed integer
    using i16  = int16_t;   ///< 16-bit signed integer
    using i32  = int32_t;   ///< 32-bit signed integer
    using i64  = int64_t;   ///< 64-bit signed integer
    using ui8  = uint8_t;   ///< 8-bit unsigned integer
    using ui16 = uint16_t;  ///< 16-bit unsigned integer
    using ui32 = uint32_t;  ///< 32-bit unsigned integer
    using ui64 = uint64_t;  ///< 64-bit unsigned integer

    // Floating-point types.
    using f32 = float;   ///< 32-bit floating point value (single)
    using f64 = double;  ///< 64-bit floating point value (double)

    // GLM Vector & Matrix
    template <size_t Length, typename Type>
    using vec = glm::vec<Length, Type, glm::qualifier::aligned_highp>;

    // GLM Vector & Matrix
    template <size_t Columns, size_t Rows, typename Type>
    using mat = glm::mat<Columns, Rows, Type, glm::qualifier::aligned_highp>;

    // Integral vector types.
    using i8v2 = vec<2, i8>;
    using i8v3 = vec<3, i8>;
    using i8v4 = vec<4, i8>;

    using i16v2 = vec<2, i16>;
    using i16v3 = vec<3, i16>;
    using i16v4 = vec<4, i16>;

    using i32v2 = vec<2, i32>;
    using i32v3 = vec<3, i32>;
    using i32v4 = vec<4, i32>;

    using i64v2 = vec<2, i64>;
    using i64v3 = vec<3, i64>;
    using i64v4 = vec<4, i64>;

    using ui8v2 = vec<2, ui8>;
    using ui8v3 = vec<3, ui8>;
    using ui8v4 = vec<4, ui8>;

    using ui16v2 = vec<2, ui16>;
    using ui16v3 = vec<3, ui16>;
    using ui16v4 = vec<4, ui16>;

    using ui32v2 = vec<2, ui32>;
    using ui32v3 = vec<3, ui32>;
    using ui32v4 = vec<4, ui32>;

    using ui64v2 = vec<2, ui64>;
    using ui64v3 = vec<3, ui64>;
    using ui64v4 = vec<4, ui64>;

    // Floating-point vector types.
    using f32v2 = vec<2, f32>;
    using f32v3 = vec<3, f32>;
    using f32v4 = vec<4, f32>;

    using f64v2 = vec<2, f64>;
    using f64v3 = vec<3, f64>;
    using f64v4 = vec<4, f64>;

    // Integral matrix types.
    using i8m2 = mat<2, 2, i8>;
    using i8m3 = mat<3, 3, i8>;
    using i8m4 = mat<4, 4, i8>;

    using i16m2 = mat<2, 2, i16>;
    using i16m3 = mat<3, 3, i16>;
    using i16m4 = mat<4, 4, i16>;

    using i32m2 = mat<2, 2, i32>;
    using i32m3 = mat<3, 3, i32>;
    using i32m4 = mat<4, 4, i32>;

    using i64m2 = mat<2, 2, i64>;
    using i64m3 = mat<3, 3, i64>;
    using i64m4 = mat<4, 4, i64>;

    using ui8m2 = mat<2, 2, ui8>;
    using ui8m3 = mat<3, 3, ui8>;
    using ui8m4 = mat<4, 4, ui8>;

    using ui16m2 = mat<2, 2, ui16>;
    using ui16m3 = mat<3, 3, ui16>;
    using ui16m4 = mat<4, 4, ui16>;

    using ui32m2 = mat<2, 2, ui32>;
    using ui32m3 = mat<3, 3, ui32>;
    using ui32m4 = mat<4, 4, ui32>;

    using ui64m2 = mat<2, 2, ui64>;
    using ui64m3 = mat<3, 3, ui64>;
    using ui64m4 = mat<4, 4, ui64>;

    // Floating-point matrix types.
    using f32m2 = mat<2, 2, f32>;
    using f32m3 = mat<3, 3, f32>;
    using f32m4 = mat<4, 4, f32>;

    using f64m2 = mat<2, 2, f64>;
    using f64m3 = mat<3, 3, f64>;
    using f64m4 = mat<4, 4, f64>;
}  // namespace nbs

#endif  // N_BODY_SIM_TYPES_HPP

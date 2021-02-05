#ifndef N_BODY_SIM_TYPES_H
#define N_BODY_SIM_TYPES_H

#pragma once

#include <cstdint>

// Decorators

// Indicates a parameter that will be initialised inside function call, including any necesasry memory allocation.
#define OUT

// Types

typedef int8_t   i8;   ///< 8-bit signed integer
typedef int16_t  i16;  ///< 16-bit signed integer
typedef int32_t  i32;  ///< 32-bit signed integer
typedef int64_t  i64;  ///< 64-bit signed integer
typedef uint8_t  ui8;  ///< 8-bit unsigned integer
typedef uint16_t ui16; ///< 16-bit unsigned integer
typedef uint32_t ui32; ///< 32-bit unsigned integer
typedef uint64_t ui64; ///< 64-bit unsigned integer

typedef float f32; ///< 32-bit floating point value (single)
typedef double f64; ///< 64-bit floating point value (double)

#endif // N_BODY_SIM_TYPES_H

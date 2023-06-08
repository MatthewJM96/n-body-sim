#ifndef N_BODY_SIM_STDAFX_H
#define N_BODY_SIM_STDAFX_H

#pragma once

// Basics
#include <cstdint>
#include <cstdlib>
#include <random>

// Generics
#include <type_traits>

// Strings
#include <cstring>
#include <string>

// Streams
#include <iostream>

// GL Maths
#define GLM_ENABLE_EXPERIMENTAL

#define GLM_FORCE_AVX2
#define GLM_FORCE_SWIZZLE
#define GLM_FORCE_SIZE_T_LENGTH
#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// Our Stuff
#include "decorators.hpp"
#include "types.hpp"

#endif  // N_BODY_SIM_STDAFX_H

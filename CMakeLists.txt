cmake_minimum_required(VERSION 3.22)

message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
            "${CMAKE_BINARY_DIR}/cmake/conan.cmake"
            EXPECTED_HASH SHA256=5cdb3042632da3efff558924eecefd580a0e786863a857ca097c3d1d43df5dcd
            TLS_VERIFY ON
            STATUS DOWNLOAD_STATUS)

list(GET DOWNLOAD_STATUS 0 DOWNLOAD_STATUS_CODE)
if(NOT ${DOWNLOAD_STATUS_CODE} EQUAL 0)
    message(FATAL_ERROR "Error downloading conan.cmake: ${DOWNLOAD_STATUS}")
else()
    message(STATUS "Download result: ${DOWNLOAD_STATUS}")
endif()

include(${CMAKE_BINARY_DIR}/cmake/conan.cmake)

project(N_Body_Simulator)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(NBS_FAST_DEBUG "Whether to compile debug builds with O1 optimisation." OFF)
option(NBS_SUPER_FAST_DEBUG "Whether to compile debug builds with O2 optimisation." OFF)
option(NBS_HYPER_FAST_DEBUG "Whether to compile debug builds with O3 optimisation." OFF)

option(NBS_ENABLE_ADDRESS_SANITIZER "Whether to compile with address sanitizer." OFF)
option(NBS_ENABLE_THREAD_SANITIZER "Whether to compile with thread sanitizer." OFF)
option(NBS_ENABLE_MEMORY_SANITIZER "Whether to compile with memory sanitizer." OFF)
option(NBS_ENABLE_GPERF_PROFILER "Whether to compile with GPerf profiler." OFF)

if (UNIX OR MINGW)
    set(CMAKE_CXX_FLAGS "-pthread -mavx2")
    if (NBS_ENABLE_ADDRESS_SANITIZER)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
        set(CMAKE_CXX_FLAGS_RELEASE "-g -O3")
        set(CMAKE_CXX_FLAGS_DEBUG "-g")
    elseif(NBS_ENABLE_THREAD_SANITIZER)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -fno-omit-frame-pointer")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=thread -fno-omit-frame-pointer")
        set(CMAKE_CXX_FLAGS_RELEASE "-g -O2")
        set(CMAKE_CXX_FLAGS_DEBUG "-g")
    elseif(NBS_ENABLE_MEMORY_SANITIZER)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory -fsanitize-memory-track-origins -fPIE -fno-omit-frame-pointer")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=memory -fsanitize-memory-track-origins -pie -fno-omit-frame-pointer")
            set(CMAKE_CXX_FLAGS_RELEASE "-g -O2")
            set(CMAKE_CXX_FLAGS_DEBUG "-g")
    elseif(NBS_ENABLE_GPERF_PROFILER)
        set(CMAKE_CXX_FLAGS_RELEASE "-g -O3")
        set(CMAKE_CXX_FLAGS_DEBUG "-g")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wold-style-cast -Wshadow -Wformat=2 -Werror")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3")
        set(CMAKE_CXX_FLAGS_DEBUG "-g")
    endif()

    if (NBS_FAST_DEBUG)
        set(CMAKE_CXX_FLAGS_DEBUG "-g -O1")
    elseif (NBS_SUPER_FAST_DEBUG)
        set(CMAKE_CXX_FLAGS_DEBUG "-g -O2")
    elseif (NBS_HYPER_FAST_DEBUG)
        set(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
    endif()

    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wduplicated-cond -Wuseless-cast")
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmodules-ts -fconcepts-diagnostics-depth=3")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wused-but-marked-unused")
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmodules -fprebuilt-module-path=${PROJECT_BINARY_DIR}/clm/${CMAKE_BUILD_TYPE}/")
    endif()
elseif (MSVC)
    set(CMAKE_CXX_FLAGS "/W3 /WX /EHsc /Zc:preprocessor /arch:AVX2")
    set(CMAKE_CXX_FLAGS_RELEASE "/MD /O2")
    set(CMAKE_CXX_FLAGS_DEBUG "/MDd /Od /Zi")

    # Don't error on secure warnings.
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS=1)
endif()

# if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmodules-ts")
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmodules -fprebuilt-module-path=${PROJECT_BINARY_DIR}/clm/${CMAKE_BUILD_TYPE}/")
# endif()

set(NBS_Requirements
    glm/0.9.9.8
    glew/2.2.0
    glfw/3.3.8
    imgui/1.89.4
    implot/0.14
)


conan_cmake_configure(REQUIRES ${NBS_Requirements}
                      GENERATORS cmake_find_package)

conan_cmake_autodetect(settings)

conan_cmake_install(PATH_OR_REFERENCE .
                    BUILD missing
                    REMOTE conancenter
                    SETTINGS ${settings})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}")

find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(imgui REQUIRED)
find_package(implot REQUIRED)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(DEBUG=1)
endif()

add_executable(NBS
    "${PROJECT_SOURCE_DIR}/src/stdafx.cpp"
    "${PROJECT_SOURCE_DIR}/tests/backends/imgui_impl_glfw.cpp"
    "${PROJECT_SOURCE_DIR}/tests/backends/imgui_impl_opengl3.cpp"
    "${PROJECT_SOURCE_DIR}/tests/main.cpp"
)

target_precompile_headers(NBS
    PUBLIC
        include/stdafx.h
)

set(NBS_Include_Dirs
    ${glm_INCLUDE_DIRS}
    ${glfw_INCLUDE_DIRS}
    ${imgui_INCLUDE_DIRS}
    ${implot_INCLUDE_DIRS}
)

# TODO(Matthew): remove the "${PROJECT_SOURCE_DIR}" line and properly build
#                engine and tests separately.
target_include_directories(NBS
    PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}"
)

target_include_directories(NBS
    SYSTEM
    PUBLIC
    ${NBS_Include_Dirs}
)

# if (UNIX OR MINGW)
#     if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#         if(CMAKE_BUILD_TYPE STREQUAL "Release" OR NBS_FAST_DEBUG)
#             target_include_directories(NBS
#                 SYSTEM
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/GCC-x64-Release/include"
#             )

#             target_link_directories(NBS
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/GCC-x64-Release/lib"
#             )
#         else()
#             target_include_directories(NBS
#                 SYSTEM
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/GCC-x64-Debug/include"
#             )

#             target_link_directories(NBS
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/GCC-x64-Debug/lib"
#             )
#         endif()
#     elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#         if(CMAKE_BUILD_TYPE STREQUAL "Release" OR NBS_FAST_DEBUG)
#             target_include_directories(NBS
#                 SYSTEM
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/Clang-x64-Release/include"
#             )

#             target_link_directories(NBS
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/Clang-x64-Release/lib"
#             )
#         else()
#             target_include_directories(NBS
#                 SYSTEM
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/Clang-x64-Debug/include"
#             )

#             target_link_directories(NBS
#                 PUBLIC
#                 "${PROJECT_SOURCE_DIR}/deps/Clang-x64-Debug/lib"
#             )
#         endif()
#     endif()
# else()
#     if(CMAKE_BUILD_TYPE STREQUAL "Release" OR NBS_FAST_DEBUG)
#         target_include_directories(NBS
#             SYSTEM
#             PUBLIC
#             "${PROJECT_SOURCE_DIR}/deps/Release/include"
#         )

#         target_link_directories(NBS
#             PUBLIC
#             "${PROJECT_SOURCE_DIR}/deps/Release/lib"
#         )
#     else()
#         target_include_directories(NBS
#             SYSTEM
#             PUBLIC
#             "${PROJECT_SOURCE_DIR}/deps/Debug/include"
#         )

#         target_link_directories(NBS
#             PUBLIC
#             "${PROJECT_SOURCE_DIR}/deps/Debug/lib"
#         )
#     endif()
# endif()

target_link_libraries(NBS GLEW::GLEW)
if (UNIX OR MINGW)
    target_link_libraries(NBS GL)
elseif (MSVC)
    target_link_libraries(NBS OpenGL32)
endif()

target_link_libraries(NBS glfw::glfw imgui::imgui implot::implot)

if(NBS_ENABLE_GPERF_PROFILER)
    target_link_libraries(NBS profiler)
endif()

// Memory Management
#define CALLER_DELETE
#define CALLEE_DELETE

// Data Direction
#define OUT
#define IN

#if defined(DEBUG)
#  define debug_printf(...) printf(__VA_ARGS__)
#else
#  define debug_printf(...)
#  if !defined(NDEBUG)
#    define NDEBUG
#  endif  // !defined(NDEBUG)
#endif

#define NBS_NON_COPYABLE(TYPE)                                                         \
  TYPE(const TYPE& rhs)            = delete;                                           \
  TYPE& operator=(const TYPE& rhs) = delete

#define NBS_COPYABLE(TYPE)                                                             \
  TYPE(const TYPE& rhs) { *this = rhs; }                                               \
  TYPE& operator=(const TYPE& rhs)

#define NBS_NON_MOVABLE(TYPE)                                                          \
  TYPE(TYPE&& rhs)            = delete;                                                \
  TYPE& operator=(TYPE&& rhs) = delete

#define NBS_MOVABLE(TYPE)                                                              \
  TYPE(TYPE&& rhs) { *this = std::forward<TYPE>(rhs); }                                \
  TYPE& operator=(TYPE&& rhs)

#if defined(NBS_COMPILER_GCC) || defined(NBS_COMPILER_CLANG)
#  define NBS_PACKED_STRUCT(DECL) DECL __attribute__((packed))
#else  // defined(NBS_COMPILER_GCC) || defined(NBS_COMPILER_CLANG)
#  define NBS_PACKED_STRUCT(DECL) __pragma(pack(push, 1)) DECL __pragma(pack(pop))
#endif

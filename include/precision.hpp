#if !defined(NBS_PRECISION)
#  if defined(NBS_USE_DOUBLE_PRECISION)
#    define NBS_PRECISION f64
#  else
#    define NBS_PRECISION f32
#  endif
#endif

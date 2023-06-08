#if !defined(NBS_PRECISION)
#  if defined(NBS_USE_DOUBLE_PRECISION)
#    define NBS_PRECISION nbs::f64
#  else
#    define NBS_PRECISION nbs::f32
#  endif
#endif

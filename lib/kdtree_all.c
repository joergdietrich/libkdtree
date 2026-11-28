/* The installed library compiles all sources as one translation unit.
   The functions in kdtree_internal.h become static, so neither the
   static nor the shared library exports anything beyond kdtree.h. */

#define KD_INTERNAL static

#include "pqueue.c"
#include "kdtree_common.c"
#include "kdtree_cartesian.c"
#include "kdtree_spherical.c"

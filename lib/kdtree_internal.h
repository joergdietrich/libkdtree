/* Functions shared between the library's source files and used by
   its tests. Not installed; programs use kdtree.h. */

#ifndef _KDTREE_INTERNAL_H_
#define _KDTREE_INTERNAL_H_

#include <pthread.h>
#include <string.h>

#include "kdtree.h"

/* kdtree_all.c defines KD_INTERNAL as static, so these functions stay
   out of the installed library's symbol table. */
#ifndef KD_INTERNAL
#define KD_INTERNAL
#endif

/* Utilities */
KD_INTERNAL void *kd_malloc(size_t size, const char *msg);
KD_INTERNAL int kd_isleaf(struct kdNode *n);
KD_INTERNAL float kd_sqr(float a);
KD_INTERNAL float kd_min(float x, float y);

/* Spherical rectangle tests */
KD_INTERNAL int kd_sph_isPointInRect(struct kdNode *node, float *min,
                                     float *max);
KD_INTERNAL int kd_sph_isRectInRect(struct kdNode *node, float *min,
                                    float *max);

/* Result queues */
KD_INTERNAL struct pqueue *kd_pqinit(struct pqueue *q, uint32_t n);
KD_INTERNAL int kd_pqinsert(struct pqueue *q, struct resItem *d);
KD_INTERNAL int kd_insertResTree(struct kdNode *node, struct pqueue *res);
KD_INTERNAL struct pqueue *kd_qinit(unsigned int *q);

#endif                          /* _KDTREE_INTERNAL_H_ */

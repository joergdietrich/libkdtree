/* Functions shared between the library's source files and used by
   its tests. Not installed; programs use kdtree.h. */

#ifndef _KDTREE_INTERNAL_H_
#define _KDTREE_INTERNAL_H_

#include <pthread.h>
#include <string.h>

#include "kdtree.h"

/* Utilities */
void *kd_malloc(size_t size, const char *msg);
int kd_isleaf(struct kdNode *n);
float kd_sqr(float a);
float kd_min(float x, float y);

/* A single node outside any tree, freed with kd_destroyTree() */
struct kdNode *kd_allocNode(struct kd_point *points, unsigned long pivot,
                            float *min, float *max, int axis, int dim);

/* Spherical rectangle tests */
int kd_sph_isPointInRect(struct kdNode *node, float *min, float *max);
int kd_sph_isRectInRect(struct kdNode *node, float *min, float *max);

/* Result queues */
struct pqueue *pqinit(struct pqueue *q, uint32_t n);
int pqinsert(struct pqueue *q, struct resItem *d);
int kd_insertResTree(struct kdNode *node, struct pqueue *res);

#endif                          /* _KDTREE_INTERNAL_H_ */

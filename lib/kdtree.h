/*! \file kdtree.h
 * \brief Include file for the kdtree library 
 */

#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>


typedef struct kd_point {
    float *point;
    void *data;
} kd_point;

/*!
 * \struct kdNode 
 * \brief kd-tree node structure definition 
 */
typedef struct kdNode {
    float *location;            /*!<vector to the node's location */
    float *min;                 /*!<vector to the min coordinates of the
                                 * hyperrectangle */
    float *max;                 /*!<vector to the max coordinates of the
                                 * hyperrectangle */
    int split;                  /*!<axis along which the tree bifurcates */
    void *data;                 /*!<pointer to an optional data structure */
    struct kdNode *left;        /*!<the left child of the tree node */
    struct kdNode *right;       /*!<the right child of the tree node */
} kdNode;

/*!
 * \struct resItem
 * \brief result items, member of a priority queue 
 */
typedef struct resItem {
    struct kdNode *node;        /*!<pointer to a kdNode */
    float dist_sq;              /*!<distance squared as the priority */
} resItem;

/*!
 * \struct pqueue
 * \brief priority queue (min-max heap)
 */
typedef struct pqueue {
    uint32_t size;              /*!<current length of the queue */
    uint32_t avail;             /*!<currently allocated queue elements */
    uint32_t step;              /*!<step size in which new elements are allocated */
    struct resItem **d;         /*!<pointer to an array of result items */
} pqueue;

/*!
 * \struct kd_thread_data
 * \brief arguments passed to the threaded kd-Tree construction
 */
typedef struct kd_thread_data {
    int max_threads;
    struct kd_point *points;
    unsigned long nPoints;
    void *(*constr) (void *);
    void (*destr) (void *);
    float *min;
    float *max;
    int depth;
    int dim;
} kd_thread_data;


#define KD_ORDERED (1)
#define KD_UNORDERED (0)

/* functions for the priority queue */
struct pqueue *pqinit(struct pqueue *q, uint32_t n);
int pqinsert(struct pqueue *q, struct resItem *d);
struct resItem **pqremove_min(struct pqueue *q, struct resItem **d);
struct resItem **pqremove_max(struct pqueue *q, struct resItem **d);
struct resItem **pqpeek_min(struct pqueue *q, struct resItem **d);
struct resItem **pqpeek_max(struct pqueue *q, struct resItem **d);


/* general utility functions */
void *kd_malloc(size_t size, const char *msg);
int kd_isleaf(struct kdNode *n);
/* Cartesian */
float kd_dist_sq(float *x, float *y, int dim);
float kd_sqr(float a);
float kd_min(float x, float y);
/* Spherical */
float kd_sph_orth_dist(float *p1, float *p2, int split);
float kd_sph_dist_sq(float *x, float *y);
float kd_sph_dist(float *x, float *y);
float kd_sph_bearing(float *p1, float *p2);
float kd_sph_xtd(float *p1, float *p2, float *p3);

/* helper functions for debugging */
void kd_printNode(struct kdNode *node, void (*printData) (void *));
void kd_printTree(struct kdNode *node);

/* Functions for building and destroying trees */
void kd_freeNode(kdNode * node, void (*destr) (void *));
struct kdNode *kd_allocNode(struct kd_point *points, unsigned long pivot,
                            float *min, float *max, int dim, int axis);
void kd_destroyTree(struct kdNode *node, void (*destr) (void *));
struct kd_thread_data *kd_buildArg(struct kd_point *points,
                                   unsigned long nPoints,
                                   void *(*constr) (void *),
                                   void (*destr) (void *), float *min,
                                   float *max, int depth, int max_threads,
                                   int dim);
struct kdNode *kd_buildTree(struct kd_point *points, unsigned long nPoints,
                            void *(*constr) (void *), void (*destr) (void *),
                            float *min, float *max, int dim, int max_threads);
void *kd_doBuildTree(void *threadarg);
struct kdNode *kd_sph_buildTree(struct kd_point *points,
                                unsigned long nPoints,
                                void *(*constr) (void *),
                                void (*destr) (void *), float *min,
                                float *max, int max_threads);
void *kd_sph_doBuildTree(void *threadarg);

/* Functions for range searches 
 * Cartesian
 */
int kd_isPointInRect(struct kdNode *node, float *min, float *max, int dim);
int kd_isRectInRect(struct kdNode *node, float *min, float *max, int dim);
int kd_rectOverlapsRect(struct kdNode *node, float *min, float *max, int dim);
struct pqueue *kd_ortRangeSearch(struct kdNode *node, float *min, float *max,
                                 int dim);
int kd_doOrtRangeSearch(struct kdNode *node, float *min, float *max, int dim,
                        struct pqueue *res);
struct kdNode *kd_nearest(struct kdNode *node, float *p, float *max_dist_sq,
                          int dim);
struct pqueue *kd_qnearest(struct kdNode *node, float *p,
                           float *max_dist_sq, unsigned int q, int dim);
int kd_doQnearest(struct kdNode *node, float *p,
                  float *max_dist_sq, unsigned int q, int dim,
                  struct pqueue *res);
struct pqueue *kd_range(struct kdNode *node, float *p, float *max_dist_sq,
                        int dim, int ordered);
int kd_doRange(struct kdNode *node, float *p, float *max_dist_sq,
               int dim, struct pqueue *res, int ordered);
/* spherical */
int kd_sph_isPointInRect(struct kdNode *node, float *min, float *max);
int kd_sph_isRectInRect(struct kdNode *node, float *min, float *max);
int kd_sph_rectOverlapsRect(struct kdNode *node, float *min, float *max);
struct pqueue *kd_sph_ortRangeSearch(struct kdNode *node, float *min,
                                     float *max);
int kd_sph_doOrtRangeSearch(struct kdNode *node, float *min, float *max,
                            struct pqueue *res);
struct kdNode *kd_sph_nearest(struct kdNode *node, float *p,
                              float *max_dist_sq);
struct pqueue *kd_sph_qnearest(struct kdNode *node, float *p,
                               float *max_dist_sq, unsigned int q);
int kd_sph_doQnearest(struct kdNode *node, float *p,
                      float *max_dist_sq, unsigned int q, struct pqueue *res);
struct pqueue *kd_sph_range(struct kdNode *node, float *p, float *max_dist_sq,
                            int ordered);
int kd_sph_doRange(struct kdNode *node, float *p, float *max_dist_sq,
                   struct pqueue *res, int ordered);

/* Functions for results heaps */
int kd_insertResTree(struct kdNode *node, struct pqueue *res);


#endif                          /* _KDTREE_H_ */

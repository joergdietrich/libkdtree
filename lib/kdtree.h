/*! \file kdtree.h
 * \brief Include file for the kdtree library
 */

#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>


/*!
 * \struct kd_point
 * \brief input point for building a tree
 *
 * The tree copies the coordinates, and copies the data only through
 * the constructor passed to kd_buildTree() or kd_sph_buildTree(). The
 * array of kd_points can be freed once the tree is built.
 */
typedef struct kd_point {
    float *point;               /*!<vector of the point's coordinates */
    void *data;                 /*!<optional data, handed to the constructor */
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


/*! Range searches return their results ordered by distance */
#define KD_ORDERED (1)
/*! Range searches return their results in no particular order */
#define KD_UNORDERED (0)


/* Building and destroying trees */
struct kdNode *kd_buildTree(struct kd_point *points, unsigned long nPoints,
                            void *(*constr) (void *), void (*destr) (void *),
                            float *min, float *max, int dim, int max_threads);
struct kdNode *kd_sph_buildTree(struct kd_point *points,
                                unsigned long nPoints,
                                void *(*constr) (void *),
                                void (*destr) (void *), float *min,
                                float *max, int max_threads);
void kd_destroyTree(struct kdNode *node, void (*destr) (void *));

/* Searches in Cartesian trees */
struct kdNode *kd_nearest(struct kdNode *node, float *p, float *max_dist_sq,
                          int dim);
struct pqueue *kd_qnearest(struct kdNode *node, float *p,
                           float *max_dist_sq, unsigned int q, int dim);
struct pqueue *kd_range(struct kdNode *node, float *p, float *max_dist_sq,
                        int dim, int ordered);
struct pqueue *kd_ortRangeSearch(struct kdNode *node, float *min, float *max,
                                 int dim);

/* Searches in spherical trees */
struct kdNode *kd_sph_nearest(struct kdNode *node, float *p,
                              float *max_dist_sq);
struct pqueue *kd_sph_qnearest(struct kdNode *node, float *p,
                               float *max_dist_sq, unsigned int q);
struct pqueue *kd_sph_range(struct kdNode *node, float *p, float *max_dist_sq,
                            int ordered);
struct pqueue *kd_sph_ortRangeSearch(struct kdNode *node, float *min,
                                     float *max);

/* Distances */
float kd_dist_sq(float *x, float *y, int dim);
float kd_sph_dist(float *x, float *y);
float kd_sph_dist_sq(float *x, float *y);

/* Result queues */
struct resItem **pqremove_min(struct pqueue *q, struct resItem **d);
struct resItem **pqremove_max(struct pqueue *q, struct resItem **d);
struct resItem **pqpeek_min(struct pqueue *q, struct resItem **d);
struct resItem **pqpeek_max(struct pqueue *q, struct resItem **d);
void pqfree(struct pqueue *q);

/* Debugging */
void kd_printNode(struct kdNode *node, void (*printData) (void *));
void kd_printTree(struct kdNode *node);


#endif                          /* _KDTREE_H_ */

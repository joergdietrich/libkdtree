/*!\file kdtree_cartesian.c
 * \brief Routines specific to the n-dimensional cartesian kd-tree
 *
 */

#include "kdtree_internal.h"

static int kd_doOrtRangeSearch(struct kdNode *node, float *min, float *max,
                               int dim, struct pqueue *res);
static int kd_doQnearest(struct kdNode *node, float *p, float *max_dist_sq,
                         unsigned int q, int dim, struct pqueue *res);
static int kd_doRange(struct kdNode *node, float *p, float *max_dist_sq,
                      int dim, struct pqueue *res, int ordered);

/* ********************************************************************

   general utility functions 

   ********************************************************************* */

/*!
 * \brief squared Euclidean distance between two points
 *
 * \param x the coordinates of the first point.
 *
 * \param y the coordinates of the second point.
 *
 * \param dim the number of dimensions.
 *
 * \return the squared distance, or -1 if x or y is NULL.
 */
inline float
kd_dist_sq(float *x, float *y, int dim)
{
    int i;
    float dsq = 0;

    if (!x || !y)
        return -1;

    for(i = 0; i < dim; i++)
        dsq += kd_sqr(x[i] - y[i]);
    return dsq;
}

inline float
kd_min(float x, float y)
{
    return x < y ? x : y;
}


/* *******************************************************************

   Functions for range searches
  
   *******************************************************************
*/

/* Returns 1 if node is a point in the hyperrectangle defined by
   minimum and maximum vectors min and max. */
static int
kd_isPointInRect(struct kdNode *node, float *min, float *max, int dim)
{
    int i;

    if (node == NULL)
        return 0;

    for(i = 0; i < dim; i++) {
        if (node->location[i] < min[i] || node->location[i] > max[i])
            return 0;
    }
    return 1;
}

/* Returns 1 if the hyperrectangle of node is fully contained within
   the HR described by the minimum and maximum vectors min and
   max. Returns 0 otherwise. */
static int
kd_isRectInRect(struct kdNode *node, float *min, float *max, int dim)
{
    int i;

    if (node == NULL)
        return 0;

    for(i = 0; i < dim; i++) {
        if (node->min[i] < min[i] || node->max[i] > max[i])
            return 0;
    }
    return 1;
}

/* Returns 1 if the hyperrectangle of node overlaps the HR described
   the minimum and maximum vectors min and max. Returns 0
   otherwise. */
static int
kd_rectOverlapsRect(struct kdNode *node, float *min, float *max, int dim)
{
    int i;

    if (node == NULL)
        return 0;

    for(i = 0; i < dim; i++) {
        if ((node->min[i] > max[i] || node->max[i] < min[i]))
            return 0;
    }
    return 1;
}

/*!
 * \brief Perform orthogonal range search (get all points in a
 * hyperrectangle).
 *
 * \param node the root node of tree to be searched.
 *
 * \param min a vector with the minimum positions of the corners of the
 * hyperrectangle containing the data.
 *
 * \param max a vector with the maximum positions of the corners of
 * the hyperrectangle containing the data.
 *
 * \param dim the dimension of the data.
 *
 * \return  Pointer to a priority queue, NULL in case of problems.
*/
struct pqueue *
kd_ortRangeSearch(struct kdNode *node, float *min, float *max, int dim)
{
    struct pqueue *res;

    if ((res = kd_pqinit(NULL, 1)) == NULL)
        return NULL;
    if (!kd_doOrtRangeSearch(node, min, max, dim, res)) {
        pqfree(res);
        return NULL;
    }
    return res;
}

/* This is the orthogonal range search. Returns 1 if okay, 0 in case
   of problems. */
static int
kd_doOrtRangeSearch(struct kdNode *node, float *min, float *max,
                    int dim, struct pqueue *res)
{

    if (node == NULL)
        return 1;

    if (kd_isleaf(node) && kd_isPointInRect(node, min, max, dim)) {
        return kd_insertResTree(node, res);
    } else {
        if (kd_isRectInRect(node->left, min, max, dim)) {
            if (!kd_insertResTree(node->left, res))
                return 0;
        } else {
            if (kd_rectOverlapsRect(node->left, min, max, dim))
                if (!kd_doOrtRangeSearch(node->left, min, max, dim, res))
                    return 0;
        }
        if (kd_isRectInRect(node->right, min, max, dim)) {
            if (!kd_insertResTree(node->right, res))
                return 0;
        } else {
            if (kd_rectOverlapsRect(node->right, min, max, dim))
                if (!kd_doOrtRangeSearch(node->right, min, max, dim, res))
                    return 0;
        }
    }
    return 1;
}

/*! 
 * \brief Find the nearest neighbor of a point.
 *
 * \param node the root node of the tree to be searched.
 *
 * \param p a vector to the point whose nearest neighbor is sought.
 *
 * \param max_dist_sq the square of the maximum distance to the
 * nearest neighbor.
 *
 * \param dim the dimension of the data.
 *
 * \return A pointer to node containing the nearest neighbor.
 * max_dist_sq is set to the square of the distance to the nearest
 * neigbor.
 */
struct kdNode *
kd_nearest(struct kdNode *node, float *p, float *max_dist_sq, int dim)
{
    struct kdNode *nearer, *further, *nearest, *tmp, *tmp_nearest;
    float dist_sq, tmp_dist_sq, dx;

    if (!node)
        return NULL;

    if (p[node->split] < node->location[node->split]) {
        nearer = node->left;
        further = node->right;
    } else {
        nearer = node->right;
        further = node->left;
    }
    tmp = kd_nearest(nearer, p, max_dist_sq, dim);
    if (tmp)
        nearest = tmp;
    else
        nearest = node;

    dist_sq = kd_dist_sq(nearest->location, p, dim);
    if (*max_dist_sq > dist_sq)
	*max_dist_sq = dist_sq;

    if (!further)
        return nearest;

    dx = kd_min(fabs(p[node->split] - further->min[node->split]),
                fabs(p[node->split] - further->max[node->split]));
    if (*max_dist_sq > kd_sqr(dx)) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        tmp = kd_nearest(further, p, max_dist_sq, dim);
        if (tmp)
            tmp_nearest = tmp;
        else
            tmp_nearest = further;

        tmp_dist_sq = kd_dist_sq(tmp_nearest->location, p, dim);
        if (tmp_dist_sq < dist_sq) {
            nearest = tmp_nearest;
            dist_sq = tmp_dist_sq;
            *max_dist_sq = kd_min(dist_sq, *max_dist_sq);
        }
    }
    return nearest;
}

/*! 
 * \brief Return the q nearest-neighbors to a point.
 *
 * \param node the root node of the tree to be searched.
 *
 * \param p a vector to the point whose nearest neighbors are sought.
 *
 * \param max_dist_sq the square of the maximum distance to the
 * nearest neighbors.
 *
 * \param q the maximum number of points to be returned.
 *
 * \param dim the dimension of the data.
 *
 * \return A pointer to a priority queue of the points found, or NULL
 * in case of problems.
 */
struct pqueue *
kd_qnearest(struct kdNode *node, float *p,
            float *max_dist_sq, unsigned int q, int dim)
{
    struct pqueue *res;

    if ((res = kd_qinit(&q)) == NULL)
        return NULL;
    if (!kd_doQnearest(node, p, max_dist_sq, q + 1, dim, res)) {
        pqfree(res);
        return NULL;
    }
    return res;
}

/* *
 * This is the q nearest-neighbor search.
 *
 * This is a modification of the range search in which the maximum
 * search radius is decreased to the maximum of the queue as soon as
 * the queue is filled.
 *
 * return 1 if okay, zero in case of problems
 */
static int
kd_doQnearest(struct kdNode *node, float *p,
              float *max_dist_sq, unsigned int q, int dim, struct pqueue *res)
{
    struct kdNode *nearer, *further;
    struct resItem *point, *item;
    float dist_sq, dx;

    if (!node)
        return 1;

    dist_sq = kd_dist_sq(node->location, p, dim);
    if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
        if ((point = kd_malloc(sizeof(struct resItem), "kd_doQnearest: "))
            == NULL)
            return 0;
        point->node = node;
        point->dist_sq = dist_sq;
        if (!kd_pqinsert(res, point)) {
            free(point);
            return 0;
        }
    }
    if (res->size > q) {
        pqremove_max(res, &item);
        free(item);
        if (res->size > 1) {
            /*
             * Only inspect the queue if there are items left 
             */
            pqpeek_max(res, &item);
            *max_dist_sq = item->dist_sq;
        } else {
            /*
             * Nothing was found within the max search radius 
             */
            *max_dist_sq = 0;
        }
    }

    if (p[node->split] < node->location[node->split]) {
        nearer = node->left;
        further = node->right;
    } else {
        nearer = node->right;
        further = node->left;
    }
    if (!kd_doQnearest(nearer, p, max_dist_sq, q, dim, res))
        return 0;

    if (!further)
        return 1;

    dx = kd_min(fabs(p[node->split] - further->min[node->split]),
                fabs(p[node->split] - further->max[node->split]));
    if (*max_dist_sq > kd_sqr(dx)) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        if (!kd_doQnearest(further, p, max_dist_sq, q, dim, res))
            return 0;
    }
    return 1;
}


/*!
 * \brief Perform a range search around a point.
 *
 * \param node the root node of the tree to be searched.
 *
 * \param p the location of the point around which the search is carried out .
 *
 * \param max_dist_sq the square of the radius of the hypersphere.
 *
 * \param dim the dimension of the data.  \param ordered determines
 * whether the result list should be ordered in increasing distance
 * (KD_ORDERED) or unordered (KD_UNORDERED).
 *
 * \return A pointer to a priority queue containing the points found,
 * NULL in case of problems.
 */
struct pqueue *
kd_range(struct kdNode *node, float *p, float *max_dist_sq,
         int dim, int ordered)
{
    struct pqueue *res;

    if ((res = kd_pqinit(NULL, 1)) == NULL)
        return NULL;
    if (!kd_doRange(node, p, max_dist_sq, dim, res, ordered)) {
        pqfree(res);
        return NULL;
    }
    return res;
}


/* This is the range search. Returns 1 if okay, 0 in case of problems */
static int
kd_doRange(struct kdNode *node, float *p, float *max_dist_sq,
           int dim, struct pqueue *res, int ordered)
{

    struct kdNode *nearer, *further;
    struct resItem *point;
    float dist_sq, dx;

    if (!node)
        return 1;

    dist_sq = kd_dist_sq(node->location, p, dim);
    if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
        if ((point = kd_malloc(sizeof(struct resItem), "kd_doRange:"))
            == NULL)
            return 0;
        point->node = node;
        point->dist_sq = ordered ? dist_sq : -1;
        if (!kd_pqinsert(res, point)) {
            free(point);
            return 0;
        }
    }

    if (p[node->split] < node->location[node->split]) {
        nearer = node->left;
        further = node->right;
    } else {
        nearer = node->right;
        further = node->left;
    }
    if (!kd_doRange(nearer, p, max_dist_sq, dim, res, ordered))
        return 0;

    if (!further)
        return 1;

    dx = kd_min(fabs(p[node->split] - further->min[node->split]),
                fabs(p[node->split] - further->max[node->split]));
    if (*max_dist_sq > kd_sqr(dx)) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        if (!kd_doRange(further, p, max_dist_sq, dim, res, ordered))
            return 0;
    }
    return 1;
}

/* End range searching functions */

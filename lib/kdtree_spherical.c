/*!\file kdtree_spherical.c
 * \brief Routines specific to the 2-dimensional spherical kd-tree
 *
 */

#include "kdtree_internal.h"

static int kd_sph_doOrtRangeSearch(struct kdNode *node, float *min,
                                   float *max, struct pqueue *res);
static int kd_sph_doQnearest(struct kdNode *node, float *p,
                             float *max_dist_sq, unsigned int q,
                             struct pqueue *res);
static int kd_sph_doRange(struct kdNode *node, float *p, float *max_dist_sq,
                          struct pqueue *res, int ordered);

/* ********************************************************************

   general utility functions 

   ********************************************************************* */

/* Great-circle distance by the haversine formula, which stays
   accurate for small separations. */
inline float
kd_sph_dist(float *x, float *y)
{
    double sdec, sra, h;

    if (!x || !y)
        return -1;

    sdec = sin((x[1] - y[1]) / 2);
    sra = sin((x[0] - y[0]) / 2);
    h = sdec * sdec + cos(x[1]) * cos(y[1]) * sra * sra;

    /* Rounding can push h past 1 for antipodal points */
    return 2 * asin(sqrt(fmin(h, 1.0)));
}

inline float
kd_sph_dist_sq(float *x, float *y)
{
    if (!x || !y)
        return -1;

    return kd_sqr(kd_sph_dist(x, y));
}

static float
kd_sph_bearing(float *p1, float *p2)
{
    float x, y;
    
    x = cos(p1[1]) * sin(p2[1]) - sin(p1[1]) * cos(p2[1]) * cos(p2[0] - p1[0]);
    y = sin(p2[0] - p1[0]) * cos(p2[1]);
    return atan2(y, x);
}


/* Cross-track distance of p3 from the great circle through p1 and p2 */
static float
kd_sph_xtd(float *p1, float *p2, float *p3)
{
    float d13, theta13, theta12;
    
    d13 = kd_sph_dist(p1, p3);
    theta13 = kd_sph_bearing(p1, p3);
    theta12 = kd_sph_bearing(p1, p2);
    return asin(sin(d13) * sin(theta13 - theta12));
}


static float
kd_sph_orth_dist(float *p1, float *p2, int split)
{
    float q[2];

    if (split == 1) {
        q[0] = p1[0];
        q[1] = p2[1];
    } else {
        q[0] = p2[0];
        q[1] = p1[1];
    }
    return kd_sph_dist(p1, q);
}

/* end utility functions */

/* ******************************************************************
   
   Functions for building and destroying trees 

   ******************************************************************** */

/*! \brief build kd-tree structure
 *
 * \param points an array of kd_points (struct with position vector
 * and data container).
 *
 * \param nPoints the length of the points array.
 *
 * \param constr a pointer to a void *constructor() function to include
 * data container in tree; optional, can be NULL
 *
 * \param destr a pointer to a void destructor() function to free()
 * data container in the tree; optional, can be NULL, but should be
 * given if the constr argument is non-NULL.
 *
 * \param min a vector with the minimum positions of the corners of the
 * hyperrectangle containing the data.
 *
 * \param max a vector with the maximum positions of the corners of
 * the hyperrectangle containing the data.
 *
 * \param max_threads the maximal number of threads spawned for
 * construction of the tree. The threads will be unbalanced if this is
 * not a power of 2.
 *
 * \return root node of the tree, NULL if points is empty or memory
 * runs out
 */
struct kdNode *
kd_sph_buildTree(struct kd_point *points, unsigned long nPoints,
                 void *(*constr) (void *), void (*destr) (void *),
                 float *min, float *max, int max_threads)
{
    return kd_buildTree(points, nPoints, constr, destr, min, max, 2,
                        max_threads);
}


/* end of tree construction and destruction */


/* *******************************************************************

   Functions for range searches 
  
   *******************************************************************
*/

/* A query rectangle with min[0] < 0 wraps around RA = 0 and sees
   RA > pi as negative. Returns ra in the frame of that rectangle. */
static double
_sph_queryRA(double ra, float *min)
{
    return (min[0] < 0 && ra > M_PI) ? ra - 2 * M_PI : ra;
}

/* Maps the RA extent of node into the frame of the query rectangle
   with lower corner min. An extent that straddles the wrap at RA = pi
   splits in two. Writes the intervals to lo and hi and returns their
   number. */
static int
_sph_nodeRA(struct kdNode *node, float *min, double *lo, double *hi)
{
    if (min[0] < 0 && node->min[0] <= M_PI && node->max[0] > M_PI) {
        lo[0] = node->min[0];
        hi[0] = M_PI;
        lo[1] = -M_PI;
        hi[1] = node->max[0] - 2 * M_PI;
        return 2;
    }
    lo[0] = _sph_queryRA(node->min[0], min);
    hi[0] = _sph_queryRA(node->max[0], min);
    return 1;
}

/* Returns 1 if node is a point in the hyperrectangle defined by
   minimum and maximum vectors min and max. */
int
kd_sph_isPointInRect(struct kdNode *node, float *min, float *max)
{
    double ra;

    if (node == NULL)
        return 0;

    if (node->location[1] < min[1] || node->location[1] > max[1])
        return 0;

    ra = _sph_queryRA(node->location[0], min);
    return ra >= min[0] && ra <= max[0];
}

/* Returns 1 if the hyperrectangle of node is fully contained within
   the HR described by the minimum and maximum vectors min and
   max. Returns 0 otherwise. */
int
kd_sph_isRectInRect(struct kdNode *node, float *min, float *max)
{
    double lo[2], hi[2];
    int i, n;

    if (node == NULL)
        return 0;

    if (node->min[1] < min[1] || node->max[1] > max[1])
        return 0;

    n = _sph_nodeRA(node, min, lo, hi);
    for(i = 0; i < n; i++) {
        if (lo[i] < min[0] || hi[i] > max[0])
            return 0;
    }
    return 1;
}

/* Returns 1 if the hyperrectangle of node overlaps the HR described
   by the minimum and maximum vectors min and max. Returns 0
   otherwise. */
static int
kd_sph_rectOverlapsRect(struct kdNode *node, float *min, float *max)
{
    double lo[2], hi[2];
    int i, n;

    if (node == NULL)
        return 0;

    if (node->min[1] > max[1] || node->max[1] < min[1])
        return 0;

    n = _sph_nodeRA(node, min, lo, hi);
    for(i = 0; i < n; i++) {
        if (lo[i] <= max[0] && hi[i] >= min[0])
            return 1;
    }
    return 0;
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
 * \return  Pointer to a priority queue, NULL in case of problems.
 *
 * A rectangle that crosses RA = 0 is given with min[0] < 0. RA values
 * above pi then count as negative.
*/
struct pqueue *
kd_sph_ortRangeSearch(struct kdNode *node, float *min, float *max)
{
    struct pqueue *res;

    if ((res = pqinit(NULL, 1)) == NULL)
        return NULL;
    if (!kd_sph_doOrtRangeSearch(node, min, max, res)) {
        pqfree(res);
        return NULL;
    }
    return res;
}


/* This is the orthogonal range search. Returns 1 if okay, 0 in case
   of problems. */
static int
kd_sph_doOrtRangeSearch(struct kdNode *node, float *min, float *max,
                        struct pqueue *res)
{

    if (node == NULL)
        return 1;

    if (kd_isleaf(node) && kd_sph_isPointInRect(node, min, max)) {
        return kd_insertResTree(node, res);
    } else {
        if (kd_sph_isRectInRect(node->left, min, max)) {
            if (!kd_insertResTree(node->left, res))
                return 0;
        } else {
            if (kd_sph_rectOverlapsRect(node->left, min, max))
                if (!kd_sph_doOrtRangeSearch(node->left, min, max, res))
                    return 0;
        }
        if (kd_sph_isRectInRect(node->right, min, max)) {
            if (!kd_insertResTree(node->right, res))
                return 0;
        } else {
            if (kd_sph_rectOverlapsRect(node->right, min, max))
                if (!kd_sph_doOrtRangeSearch(node->right, min, max, res))
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
 * \return A pointer to node containing the nearest neighbor.
 * max_dist_sq is set to the square of the distance to the nearest
 * neigbor.
 */
struct kdNode *
kd_sph_nearest(struct kdNode *node, float *p, float *max_dist_sq)
{
    struct kdNode *nearer, *further, *nearest, *tmp, *tmp_nearest;
    float dist_sq, tmp_dist_sq, dx;
    float p1[2], p2[2];

    if (!node)
        return NULL;

    /* This is a "<" sign, not the kd_sph_lt function because the tree
       is constructed assuming a Cartesian geometry. */
    if (p[node->split] < node->location[node->split]) {
        nearer = node->left;
        further = node->right;
    } else {
        nearer = node->right;
        further = node->left;
    }
    tmp = kd_sph_nearest(nearer, p, max_dist_sq);
    if (tmp)
        nearest = tmp;
    else
        nearest = node;

    dist_sq = kd_sph_dist_sq(nearest->location, p);
    if (*max_dist_sq > dist_sq)
	*max_dist_sq = dist_sq;

    if (!further)
        return nearest;

    p1[0] = further->min[0];
    p1[1] = further->max[1];
    p2[0] = further->max[0];
    p2[1] = further->min[1];
    
    if (node->split == 0)
	/* We need to know the distance between the point p and the
	   vertical lines delineating the rectangle. These lines are
	   great circles so we need to know the cross-track
	   distance. */
	dx = kd_min(kd_sqr(kd_sph_xtd(further->min, p1, p)),
		    kd_sqr(kd_sph_xtd(p2, further->max, p)));
    else
	/* We need to know the distance between the point p and the
	   horizontal lines delineating the rectangle. Since the
	   latitudes are not great circles, the shortest distance runs
	   along the longitude of point p. */
	dx = kd_min(kd_sqr(kd_sph_orth_dist(further->min, p, 1)),
		    kd_sqr(kd_sph_orth_dist(further->max, p, 1)));

    if (*max_dist_sq > dx) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        tmp = kd_sph_nearest(further, p, max_dist_sq);
        if (tmp)
            tmp_nearest = tmp;
        else
            tmp_nearest = further;

        tmp_dist_sq = kd_sph_dist_sq(tmp_nearest->location, p);
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
 * \param q the maximum number of points to be retured.
 *
 * \return A pointer to a priority queue of the points found, or NULL
 * in case of problems.
 */
struct pqueue *
kd_sph_qnearest(struct kdNode *node, float *p,
                float *max_dist_sq, unsigned int q)
{
    struct pqueue *res;

    if ((res = pqinit(NULL, q + 2)) == NULL)
        return NULL;
    if (!kd_sph_doQnearest(node, p, max_dist_sq, q + 1, res)) {
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
kd_sph_doQnearest(struct kdNode *node, float *p, float *max_dist_sq,
                  unsigned int q, struct pqueue *res)
{
    struct kdNode *nearer, *further;
    struct resItem *point, *item;
    float dist_sq, dx;
    float p1[2], p2[2];

    if (!node)
        return 1;

    dist_sq = kd_sph_dist_sq(node->location, p);
    if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
        if ((point = kd_malloc(sizeof(struct resItem), "kd_sph_doQnearest: "))
            == NULL)
            return 0;
        point->node = node;
        point->dist_sq = dist_sq;
        if (!pqinsert(res, point)) {
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
    if (!kd_sph_doQnearest(nearer, p, max_dist_sq, q, res))
        return 0;

    if (!further)
        return 1;

    p1[0] = further->min[0];
    p1[1] = further->max[1];
    p2[0] = further->max[0];
    p2[1] = further->min[1];
    
    if (node->split == 0)
	/* We need to know the distance between the point p and the
	   vertical lines delineating the rectangle. These lines are
	   great circles so we need to know the cross-track
	   distance. */
	dx = kd_min(kd_sqr(kd_sph_xtd(further->min, p1, p)),
		    kd_sqr(kd_sph_xtd(p2, further->max, p)));
    else
	/* We need to know the distance between the point p and the
	   horizontal lines delineating the rectangle. Since the
	   latitudes are not great circles, the shortest distance runs
	   along the longitude of point p. */
	dx = kd_min(kd_sqr(kd_sph_orth_dist(further->min, p, 1)),
		    kd_sqr(kd_sph_orth_dist(further->max, p, 1)));

    if (*max_dist_sq > dx) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        if (!kd_sph_doQnearest(further, p, max_dist_sq, q, res))
            return 0;
        dist_sq = kd_sph_dist_sq(node->location, p);

        if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
            if ((point =
                 kd_malloc(sizeof(struct resItem), "kd_sph_doQnearest: "))
                == NULL)
                return 0;
            point->node = node;
            point->dist_sq = dist_sq;
            if (!pqinsert(res, point)) {
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
 * \param ordered determines whether the result list should be ordered
 * in increasing distance (KD_ORDERED) or unordered (KD_UNORDERED).
 *
 * \return A pointer to a priority queue containing the points found,
 * NULL in case of problems.
 */
struct pqueue *
kd_sph_range(struct kdNode *node, float *p, float *max_dist_sq, int ordered)
{
    struct pqueue *res;

    if ((res = pqinit(NULL, 1)) == NULL)
        return NULL;
    if (!kd_sph_doRange(node, p, max_dist_sq, res, ordered)) {
        pqfree(res);
        return NULL;
    }
    return res;
}


/* This is the range search. Returns 1 if okay, 0 in case of problems */
static int
kd_sph_doRange(struct kdNode *node, float *p, float *max_dist_sq,
               struct pqueue *res, int ordered)
{

    struct kdNode *nearer, *further;
    struct resItem *point;
    float dist_sq, dx;
    float p1[2], p2[2];

    if (!node)
        return 1;

    dist_sq = kd_sph_dist_sq(node->location, p);
    if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
        if ((point = kd_malloc(sizeof(struct resItem), "kd_sph_doRange:"))
            == NULL)
            return 0;
        point->node = node;
        point->dist_sq = ordered ? dist_sq : -1;
        if (!pqinsert(res, point)) {
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
    if (!kd_sph_doRange(nearer, p, max_dist_sq, res, ordered))
        return 0;

    if (!further)
        return 1;

    p1[0] = further->min[0];
    p1[1] = further->max[1];
    p2[0] = further->max[0];
    p2[1] = further->min[1];
    
    if (node->split == 0)
	/* We need to know the distance between the point p and the
	   vertical lines delineating the rectangle. These lines are
	   great circles so we need to know the cross-track
	   distance. */
	dx = kd_min(kd_sqr(kd_sph_xtd(further->min, p1, p)),
		    kd_sqr(kd_sph_xtd(p2, further->max, p)));
    else
	/* We need to know the distance between the point p and the
	   horizontal lines delineating the rectangle. Since the
	   latitudes are not great circles, the shortest distance runs
	   along the longitude of point p. */
	dx = kd_min(kd_sqr(kd_sph_orth_dist(further->min, p, 1)),
		    kd_sqr(kd_sph_orth_dist(further->max, p, 1)));

    if (*max_dist_sq > dx) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        if (!kd_sph_doRange(further, p, max_dist_sq, res, ordered))
            return 0;
        dist_sq = kd_sph_dist_sq(node->location, p);

        if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
            if ((point =
                 kd_malloc(sizeof(struct resItem), "kd_sph_doRange: "))
                == NULL)
                return 0;
            point->node = node;
            point->dist_sq = ordered ? dist_sq : -1;
            if (!pqinsert(res, point)) {
                free(point);
                return 0;
            }
        }

    }
    return 1;
}

/* End range searching functions */

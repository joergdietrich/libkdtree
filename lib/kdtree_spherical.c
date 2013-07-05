/*!\file kdtree_spherical.c
 * \brief Routines specific to the 2-dimensional spherical kd-tree
 *
 */

#include "kdtree.h"
#include "pqueue.h"

/* ********************************************************************

   general utility functions 

   ********************************************************************* */
inline float
kd_sph_dist_sq(float *x, float *y)
{
    float ds;
    /*float arg;*/

    if (!x || !y)
        return -1;

    ds = acos(sin(x[1]) * sin(y[1]) +
               cos(x[1]) * cos(y[1]) * cos(x[0] - y[0]));
    /*    arg = sqrt(kd_sqr(sin((x[1] - y[1])/2)) + cos(x[1]) * cos(y[1]) *
	       kd_sqr(sin((x[0] - y[0])/2)));
	       ds = 2 * asin(arg);
    */
    return kd_sqr(ds);
}

inline float
kd_sph_dist(float *x, float *y)
{
    float ds;

    if (!x || !y)
        return -1;

    ds = acos(sin(x[1]) * sin(y[1]) +
               cos(x[1]) * cos(y[1]) * cos(x[0] - y[0]));
    return ds;
}

float
kd_sph_bearing(float *p1, float *p2)
{
    float x, y;
    
    x = cos(p1[1]) * sin(p2[1]) - sin(p1[1]) * cos(p2[1]) * cos(p2[0] - p1[0]);
    y = sin(p2[0] - p1[0]) * cos(p2[1]);
    return atan2(y, x);
}


/*! \brief Compute cross-track distance
 * \param p1 a point
 *
 * \param p2 a second point defining a great circle from p1
 *
 * \param p3 a third point whose distance from the great circle we compute
 *
 * \return distance of p3 from the great circle connecting p1 and p2. 
 */
float 
kd_sph_xtd(float *p1, float *p2, float *p3)
{
    float d13, theta13, theta12;
    
    d13 = kd_sph_dist(p1, p3);
    theta13 = kd_sph_bearing(p1, p3);
    theta12 = kd_sph_bearing(p1, p2);
    return asin(sin(d13) * sin(theta13 - theta12));
}


float
kd_sph_orth_dist(float *p1, float *p2, int split)
{

    float ra2, dec2;
    float dx;

    if (split == 1) {
        ra2 = p1[0];
        dec2 = p2[1];
    } else {
        ra2 = p2[0];
        dec2 = p1[1];
    }
    dx = acos(sin(p1[1]) * sin(dec2) +
              cos(p1[1]) * cos(dec2) * cos(p1[0] - ra2));
    return dx;
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
 * \return root node of the tree
 */
struct kdNode *
kd_sph_buildTree(struct kd_point *points, unsigned long nPoints,
                 void *(*constr) (void *), void (*destr) (void *),
                 float *min, float *max, int max_threads)
{
    struct kd_thread_data *my_data;
    struct kdNode *tree;

    my_data = kd_buildArg(points, nPoints, constr, destr, min, max, 0,
                          max_threads, 2);
    tree = kd_doBuildTree(my_data);
    free(my_data);
    return tree;
}


/* end of tree construction and destruction */


/* *******************************************************************

   Functions for range searches 
  
   *******************************************************************
*/

/* Returns 1 if node is a point in the hyperrectangle defined by
   minimum and maximum vectors min and max. */
int
kd_sph_isPointInRect(struct kdNode *node, float *min, float *max)
{
    if (node == NULL)
        return 0;

    if (node->location[1] < min[1] || node->location[1] > max[1])
	return 0;

    if (min[0] < 0 && node->location[0] > M_PI) {
	if (node->location[0] - 2 * M_PI < min[0] || 
	    node->location[0] - 2 * M_PI > max[0])
	    return 0;
    } else {
	if (node->location[0] < min[0] || node->location[0] > max[0]) 
	    return 0;
    }
    
    return 1;
}

/* Returns 1 if the hyperrectangle of node is fully contained within
   the HR described by the minimum and maximum vectors min and
   max. Returns 0 otherwise. */
int
kd_sph_isRectInRect(struct kdNode *node, float *min, float *max)
{
    if (node == NULL)
        return 0;

    /* if (min[0] < 0 && node->min[0] > M_PI) {
	if (node->min[0] - 2 * M_PI < min[0] || node->max[0] > max[0])
	    return 0;
    } else {
	if (node->min[0] < min[0] || node->max[0] > max[0] )
	    return 0;
    }
    if (node->min[1] < min[1] || node->max[1] > max[1])
	return 0;
    */
    
    /* If node is fully contained all corners must be points inside
       the rectangle */ 
    struct kd_point *point;
    struct kdNode *tmpNode;
    point = kd_malloc(sizeof(kd_point), "kd_point in RectInRect");
    point->point = kd_malloc(2 * sizeof(float), 
			     "point position in RectInRect");
    point->data = NULL;
    point->point[0] = node->min[0];
    point->point[1] = node->min[1];
    tmpNode = kd_allocNode(point, 0, point->point, point->point, 0, 2);
    if (!kd_sph_isPointInRect(tmpNode, min, max)) {
	kd_destroyTree(tmpNode, NULL);	
	free(point->point);
	free(point);
	return 0;
    }
    point->point[0] = node->min[0];
    point->point[1] = node->max[1];
    tmpNode = kd_allocNode(point, 0, point->point, point->point, 0, 2);
    if (!kd_sph_isPointInRect(tmpNode, min, max)) {
	kd_destroyTree(tmpNode, NULL);	
	free(point->point);
	free(point);
	return 0;
    }
    point->point[0] = node->max[0];
    point->point[1] = node->min[1];
    tmpNode = kd_allocNode(point, 0, point->point, point->point, 0, 2);
    if (!kd_sph_isPointInRect(tmpNode, min, max)) {
	kd_destroyTree(tmpNode, NULL);	
	free(point->point);
	free(point);
	return 0;
    }
    point->point[0] = node->max[0];
    point->point[1] = node->max[1];
    tmpNode = kd_allocNode(point, 0, point->point, point->point, 0, 2);
    if (!kd_sph_isPointInRect(tmpNode, min, max)) {
	kd_destroyTree(tmpNode, NULL);	
	free(point->point);
	free(point);
	return 0;
    }
    kd_destroyTree(tmpNode, NULL);	
    free(point->point);
    free(point);

    return 1;
}

/* Returns 1 if the hyperrectangle of node overlaps the HR described
   by the minimum and maximum vectors min and max. Returns 0
   otherwise. */
int
kd_sph_rectOverlapsRect(struct kdNode *node, float *min, float *max)
{
    if (node == NULL)
        return 0;

    if (node->min[1] > max[1] || node->max[1] < min[1])
	return 0;

    if ((node->min[0]> max[0]) || (node->max[0] <  min[0]))
	return 0;

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
 * \return  Pointer to a priority queue, NULL in case of problems.
 *
 * Rectangle must not cross the meridian!
*/
struct pqueue *
kd_sph_ortRangeSearch(struct kdNode *node, float *min, float *max)
{
    struct pqueue *res;
    uint32_t i;

    if ((res = pqinit(NULL, 1)) == NULL)
        return NULL;
    if (!kd_sph_doOrtRangeSearch(node, min, max, res)) {
        for(i = 0; i < res->size; i++) {
            free(res->d[i]);
        }
        free(res->d);
        free(res);
        return NULL;
    }
    return res;
}


/* This is the orthogonal range search. Returns 1 if okay, 0 in case
   of problems. */
int
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
    uint32_t i;

    if ((res = pqinit(NULL, q + 2)) == NULL)
        return NULL;
    if (!kd_sph_doQnearest(node, p, max_dist_sq, q + 1, res)) {
        for(i = 0; i < res->size; i++) {
            free(res->d[i]);
        }
        free(res->d);
        free(res);
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
int
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
        pqinsert(res, point);
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
            pqinsert(res, point);
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
    uint32_t i;

    if ((res = pqinit(NULL, 1)) == NULL)
        return NULL;
    if (!kd_sph_doRange(node, p, max_dist_sq, res, ordered)) {
        for(i = 0; i < res->size; i++) {
            free(res->d[i]);
        }
        free(res->d);
        free(res);
        return NULL;
    }
    return res;
}


/* This is the range search. Returns 1 if okay, 0 in case of problems */
int
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
        pqinsert(res, point);
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
            pqinsert(res, point);
        }

    }
    return 1;
}

/* End range searching functions */

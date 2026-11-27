/*!\file kdtree_common.c
 * \brief Routines common to the Cartesian and spherical kd-tree versions
 *
 */

#include <limits.h>

#include "kdtree_internal.h"

/* Arguments for building a subtree, possibly in its own thread */
struct kd_thread_data {
    int max_threads;
    struct kd_point *points;    /* points of the subtree */
    unsigned long nPoints;      /* number of points */
    void *(*constr) (void *);
    struct kdNode *nodes;       /* the subtree's 2 * nPoints - 1 nodes */
    float *coords;              /* their coordinates, 3 * dim per node */
    int depth;
    int dim;
};

static void *kd_doBuildTree(void *threadarg);




/* ********************************************************************

   general utility functions 

   ********************************************************************* */

inline void *
kd_malloc(size_t size, const char *msg)
{
    void *ptr;
    if ((ptr = malloc(size)) == NULL)
        perror(msg);
    return ptr;
}

inline float
kd_sqr(float x)
{
    return x == 0.0 ? 0.0 : x * x;
}

inline int
kd_isleaf(struct kdNode *n)
{
    return (n->left || n->right) ? 0 : 1;
}

/* Returns an empty result queue for a search for the q nearest
   neighbors. The search compares against q + 1, so q is lowered by one
   if that would wrap around. Up to 1024 items are allocated up front,
   and the queue grows from there. */
struct pqueue *
kd_qinit(unsigned int *q)
{
    if (*q == UINT_MAX)
        *q = UINT_MAX - 1;
    return kd_pqinit(NULL, *q < 1024 ? *q + 2 : 1026);
}

/* end utility functions */

/* *******************************************************************
   
   Helper functions for debugging

   ******************************************************************* */

/*!
 * \brief print a node to stdout for debugging
 *
 * Prints the node's address, location, split axis, bounding box and
 * children. Only the first two coordinates are shown.
 *
 * \param node the node to print. A NULL node prints a notice to stderr.
 *
 * \param printData a function that prints the node's data, or NULL to
 * print the data pointer instead.
 */
void
kd_printNode(struct kdNode *node, void (*printData) (void *))
{
    if (!node) {
        fprintf(stderr, "Node is empty.\n");
        return;
    }
    printf("Node %p at (%f, %f)\n", (void *) node, node->location[0],
           node->location[1]);

    printf("Split axis: %d\n", node->split);
    printf("Corners: (%f, %f)\t(%f, %f)\n", node->min[0], node->min[1],
           node->max[0], node->max[1]);
    printf("Children: %p\t%p\n", (void *) node->left, (void *) node->right);
    if (printData && node->data)
        printData(node->data);
    else
        printf("Data: %p\n", node->data);

    printf("\n");

}

/*!
 * \brief print the points of a tree to stdout for debugging
 *
 * Prints the first two coordinates of each point, one point per line,
 * in tree order.
 *
 * \param node the root of the tree, or of the subtree, to print.
 */
void
kd_printTree(struct kdNode *node)
{
    if (node == NULL)
        return;

    kd_printTree(node->left);
    if (kd_isleaf(node) && node->location)
        printf("%f\t%f\n", node->location[0], node->location[1]);
    kd_printTree(node->right);
}

/* End helper functions */

/* ******************************************************************
   
   Functions for building and destroying trees 

   ******************************************************************** */
static void
_swapPoints(struct kd_point *points, unsigned long i, unsigned long j)
{
    struct kd_point tmp = points[i];

    points[i] = points[j];
    points[j] = tmp;
}

/* Median of the coordinates of three points along axis */
static float
_median3(struct kd_point *points, unsigned long i, unsigned long j,
         unsigned long k, int axis)
{
    float a = points[i].point[axis];
    float b = points[j].point[axis];
    float c = points[k].point[axis];

    if (a < b)
        return b < c ? b : (a < c ? c : a);
    return a < c ? a : (b < c ? c : b);
}

/* Restores the max-heap order along axis below root, in a heap of n
   points */
static void
_siftDown(struct kd_point *points, unsigned long root, unsigned long n,
          int axis)
{
    unsigned long child;

    while ((child = 2 * root + 1) < n) {
        if (child + 1 < n &&
            points[child].point[axis] < points[child + 1].point[axis])
            child++;
        if (!(points[root].point[axis] < points[child].point[axis]))
            return;
        _swapPoints(points, root, child);
        root = child;
    }
}

/* Sorts n points along axis in O(n log n) time for any input order */
static void
_heapSort(struct kd_point *points, unsigned long n, int axis)
{
    unsigned long i;

    for(i = n / 2; i-- > 0;)
        _siftDown(points, i, n, axis);
    for(i = n - 1; i > 0; i--) {
        _swapPoints(points, 0, i);
        _siftDown(points, 0, i, axis);
    }
}

/* Reorders the n points so that points[k] holds the point a full sort
   along axis would put there, with no greater coordinate before it and
   no smaller one after it. This is quickselect with a three-way
   partition, which keeps runs of equal coordinates linear. Like
   introselect, it sorts the remaining range once it has used 2 log2(n)
   partitions, so inputs that defeat the median-of-3 pivot cost
   O(n log n) instead of O(n^2). */
static void
_selectPoint(struct kd_point *points, unsigned long n, unsigned long k,
             int axis)
{
    unsigned long lo = 0, hi = n - 1, m;
    int rounds = 0;

    for(m = n; m > 1; m >>= 1)
        rounds += 2;

    while (lo < hi) {
        unsigned long lt = lo, i = lo, gt = hi;
        float pivot;

        if (rounds-- == 0) {
            _heapSort(points + lo, hi - lo + 1, axis);
            return;
        }
        pivot = _median3(points, lo, lo + (hi - lo) / 2, hi, axis);

        /*
         * Partition [lo, hi] into coordinates below the pivot at
         * [lo, lt), equal to it at [lt, gt] and above it at (gt, hi].
         * The pivot occurs in the range, so gt never drops below lo.
         */
        while (i <= gt) {
            float x = points[i].point[axis];

            if (x < pivot)
                _swapPoints(points, lt++, i++);
            else if (x > pivot)
                _swapPoints(points, i, gt--);
            else
                i++;
        }
        if (k < lt)
            hi = lt - 1;
        else if (k > gt)
            lo = gt + 1;
        else
            return;
    }
}

/* Returns 1 if any of the n points has a NaN coordinate */
static int
_hasNaN(struct kd_point *points, unsigned long n, int dim)
{
    unsigned long i;
    int j;

    for(i = 0; i < n; i++)
        for(j = 0; j < dim; j++)
            if (isnan(points[i].point[j]))
                return 1;
    return 0;
}

/* Returns 1 if a corner of the box from min to max has a NaN coordinate */
static int
_boxHasNaN(float *min, float *max, int dim)
{
    int j;

    for(j = 0; j < dim; j++)
        if (isnan(min[j]) || isnan(max[j]))
            return 1;
    return 0;
}

/* Points the node at its slot of 3 * dim floats, which holds its
   location, then its min corner, then its max corner. */
static void
_placeNode(struct kdNode *node, float *coords, int dim)
{
    node->location = coords;
    node->min = coords + dim;
    node->max = coords + 2 * dim;
}

/* Hands the data of every node below and including node to destr */
static void
_destroyData(struct kdNode *node, void (*destr) (void *))
{
    if (node == NULL)
        return;
    _destroyData(node->left, destr);
    _destroyData(node->right, destr);
    if (node->data)
        destr(node->data);
}

/*! \brief build kd-tree structure
 *
 * \param points an array of kd_points (struct with position vector
 * and data container). The array is reordered.
 *
 * \param nPoints the length of the points array.
 *
 * \param constr a pointer to a void *constructor() function to include
 * the data container in the tree; optional, can be NULL. If it
 * returns NULL for non-NULL data, the build fails.
 *
 * \param destr a pointer to a void destructor() function to free()
 * the data containers in the tree; optional, can be NULL, but should
 * be given if the constr argument is non-NULL.
 *
 * \param min a vector with the minimum positions of the corners of the
 * hyperrectangle containing the data.
 *
 * \param max a vector with the maximum positions of the corners of
 * the hyperrectangle containing the data.
 *
 * \param dim the dimensionality of the data.
 *
 * \param max_threads the maximal number of threads spawned for
 * construction of the tree. The threads will be unbalanced if this is
 * not a power of 2.
 *
 * \return root node of the tree, NULL if points is empty, min or max
 * is NULL, dim < 1, a point or box coordinate is NaN, constr fails or
 * memory runs out
 *
 * The tree occupies two blocks of memory, one for its 2 * nPoints - 1
 * nodes and one for their coordinates.
 */
struct kdNode *
kd_buildTree(struct kd_point *points, unsigned long nPoints,
             void *(*constr) (void *), void (*destr) (void *),
             float *min, float *max, int dim, int max_threads)
{
    struct kd_thread_data arg;
    struct kdNode *nodes;
    float *coords;
    size_t nnodes, stride;

    if (!points || nPoints == 0) {
        fprintf(stderr, "kd_buildTree: no points\n");
        return NULL;
    }
    if (!min || !max) {
        fprintf(stderr, "kd_buildTree: no bounding box\n");
        return NULL;
    }
    if (dim < 1) {
        fprintf(stderr, "kd_buildTree: dim must be at least 1\n");
        return NULL;
    }
    if (_hasNaN(points, nPoints, dim) || _boxHasNaN(min, max, dim)) {
        fprintf(stderr, "kd_buildTree: NaN coordinate\n");
        return NULL;
    }
    stride = 3 * (size_t) dim;
    if (nPoints > SIZE_MAX / (2 * stride * sizeof(float))) {
        fprintf(stderr, "kd_buildTree: too many points\n");
        return NULL;
    }
    nnodes = 2 * (size_t) nPoints - 1;

    nodes = calloc(nnodes, sizeof *nodes);
    coords = malloc(nnodes * stride * sizeof *coords);
    if (!nodes || !coords) {
        perror("kd_buildTree");
        free(nodes);
        free(coords);
        return NULL;
    }

    _placeNode(nodes, coords, dim);
    memcpy(nodes->min, min, dim * sizeof(float));
    memcpy(nodes->max, max, dim * sizeof(float));

    arg.points = points;
    arg.nPoints = nPoints;
    arg.constr = constr;
    arg.nodes = nodes;
    arg.coords = coords;
    arg.depth = 0;
    arg.max_threads = max_threads;
    arg.dim = dim;
    if (!kd_doBuildTree(&arg)) {
        fprintf(stderr, "kd_buildTree: constr failed\n");
        kd_destroyTree(nodes, destr);
        return NULL;
    }
    return nodes;
}

/* Builds a subtree into the nodes and coordinates of a->nodes and
   a->coords, laid out in pre-order. The caller has placed the root of
   the subtree and filled in its box. */
static void *
kd_doBuildTree(void *threadarg)
{
    struct kd_thread_data *a = threadarg;
    struct kd_thread_data arg[2];
    struct kdNode *node = a->nodes, *child[2];
    size_t box = a->dim * sizeof(float);
    size_t stride = 3 * (size_t) a->dim;
    pthread_t threads[2];
    int threaded[2];
    unsigned long pivot;
    int axis, i;
    void *ret;

    axis = a->depth % a->dim;
    node->split = axis;
    if (a->nPoints == 1) {
        memcpy(node->location, a->points[0].point, box);
        if (a->constr) {
            node->data = a->constr(a->points[0].data);
            if (!node->data && a->points[0].data)
                return NULL;
        }
        return node;
    }

    /*
     * The median along the split axis becomes this node, with the
     * points below it to its left and those above it to its right
     */
    pivot = a->nPoints / 2;
    _selectPoint(a->points, a->nPoints, pivot, axis);
    memcpy(node->location, a->points[pivot].point, box);

    /*
     * The left subtree of 2 * pivot - 1 nodes follows this node, and
     * the right subtree follows the left one. Each child gets this
     * node's box with one side moved to the split.
     */
    node->left = node + 1;
    node->right = node + 2 * pivot;
    for(i = 0; i < 2; i++) {
        child[i] = i ? node->right : node->left;
        arg[i] = *a;
        arg[i].nodes = child[i];
        arg[i].coords = a->coords + (child[i] - node) * stride;
        arg[i].depth = a->depth + 1;
        arg[i].max_threads = a->max_threads / 2;
        _placeNode(child[i], arg[i].coords, a->dim);
        memcpy(child[i]->min, node->min, box);
        memcpy(child[i]->max, node->max, box);
    }
    node->left->max[axis] = node->location[axis];
    node->right->min[axis] = node->location[axis];
    arg[0].nPoints = pivot;
    arg[1].points = a->points + pivot;
    arg[1].nPoints = a->nPoints - pivot;

    /*
     * Build each child in a new thread if allowed. A child whose
     * thread cannot be started is built in this thread, unless the
     * left child has already failed.
     */
    for(i = 0; i < 2; i++)
        threaded[i] = a->max_threads > 1 &&
            !pthread_create(&threads[i], NULL, kd_doBuildTree, &arg[i]);
    for(i = 0; i < 2; i++) {
        if (threaded[i]) {
            pthread_join(threads[i], &ret);
            child[i] = ret;
        } else if (i == 0 || child[0]) {
            child[i] = kd_doBuildTree(&arg[i]);
        } else {
            child[i] = NULL;
        }
    }
    return child[0] && child[1] ? node : NULL;
}

/* Allocates a single node outside any tree, with its location at
   points[pivot] and its box from min to max. Returns NULL if memory
   runs out. kd_destroyTree() frees the node. */
struct kdNode *
kd_allocNode(struct kd_point *points, unsigned long pivot,
             float *min, float *max, int axis, int dim)
{
    struct kdNode *node;
    float *coords;

    if ((node = kd_malloc(sizeof(kdNode), "kd_allocNode (node): ")) == NULL)
        return NULL;
    if ((coords = kd_malloc(3 * dim * sizeof(float),
                            "kd_allocNode (coords): ")) == NULL) {
        free(node);
        return NULL;
    }
    _placeNode(node, coords, dim);
    memcpy(node->location, points[pivot].point, dim * sizeof(float));
    memcpy(node->min, min, dim * sizeof(float));
    memcpy(node->max, max, dim * sizeof(float));
    node->left = node->right = NULL;
    node->data = NULL;
    node->split = axis;
    return node;
}

/*!
 * \brief free the kd-tree data structure
 *
 * \param node the root node of the tree to be destroyed, as returned
 * by kd_buildTree() or kd_sph_buildTree()
 *
 * \param destr optional destructor, called for the data of every node
 * that has any
 */
void
kd_destroyTree(struct kdNode *node, void (*destr) (void *))
{
    if (node == NULL)
        return;
    if (destr)
        _destroyData(node, destr);

    /* The root starts both the node and the coordinate block */
    free(node->location);
    free(node);
}

/* end of tree construction and destruction */

/* Functions dealing with result heaps */

/* Insert the sub-tree starting at node into the result heap res */
int
kd_insertResTree(struct kdNode *node, struct pqueue *res)
{
    struct resItem *point;

    if (node == NULL)
        return 1;

    if (!kd_insertResTree(node->left, res))
        return 0;
    if (kd_isleaf(node)) {
        if ((point = kd_malloc(sizeof(struct resItem), "kd_insertResTree: "))
            == NULL)
            return 0;

        point->node = node;
        point->dist_sq = -1;
        if (!kd_pqinsert(res, point)) {
            free(point);
            return 0;
        }
    }
    if (!kd_insertResTree(node->right, res))
        return 0;
    return 1;
}

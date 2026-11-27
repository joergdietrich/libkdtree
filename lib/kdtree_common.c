/*!\file kdtree_common.c
 * \brief Routines common to the Cartesian and spherical kd-tree versions
 *
 */

#include "kdtree.h"
#include "pqueue.h"


extern int pmergesort(void *base, size_t nmemb, size_t size,
                      int (*compar) (const void *, const void *, int axis),
                      int axis, int max_threads);


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

/* end utility functions */

/* *******************************************************************
   
   Helper functions for debugging

   ******************************************************************* */

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
static int
_compPoints(const void *p1, const void *p2, int axis)
{
    struct kd_point *a = (struct kd_point *) p1;
    struct kd_point *b = (struct kd_point *) p2;

    if (a->point[axis] > b->point[axis])
        return 1;
    else if (a->point[axis] == b->point[axis])
        return 0;
    else
        return -1;
}

/* Arguments for building one child of a node that splits along axis
   at split: n points starting at points, inside the parent's box with
   its upper (left child) or lower (right child) side moved to split. */
static struct kd_thread_data *
_childArg(const struct kd_thread_data *parent, struct kd_point *points,
          unsigned long n, int axis, float split, int is_right)
{
    struct kd_thread_data *arg;
    size_t box = parent->dim * sizeof(float);
    float *min = kd_malloc(box, "kd_buildTree (min): ");
    float *max = kd_malloc(box, "kd_buildTree (max): ");

    if (!min || !max) {
        free(min);
        free(max);
        return NULL;
    }
    memcpy(min, parent->min, box);
    memcpy(max, parent->max, box);
    if (is_right)
        min[axis] = split;
    else
        max[axis] = split;

    arg = kd_buildArg(points, n, parent->constr, parent->destr, min, max,
                      parent->depth + 1, parent->max_threads / 2,
                      parent->dim);
    if (!arg) {
        free(min);
        free(max);
    }
    return arg;
}

static void
_freeChildArg(struct kd_thread_data *arg)
{
    if (!arg)
        return;
    free(arg->min);
    free(arg->max);
    free(arg);
}

void *
kd_doBuildTree(void *threadarg)
{
    struct kd_thread_data *a = threadarg;
    struct kd_thread_data *arg[2];
    struct kdNode *node, *child[2];
    pthread_t threads[2];
    int threaded[2];
    unsigned long pivot;
    int axis, i;
    void *ret;

    axis = a->depth % a->dim;
    if (a->nPoints == 1) {
        if ((node = kd_allocNode(a->points, 0, a->min, a->max, axis,
                                 a->dim)) == NULL)
            return NULL;
        if (a->constr)
            node->data = a->constr(a->points[0].data);
        return node;
    }

    /*
     * Sort along the split axis, with up to max_threads threads
     */
    if (!pmergesort(a->points, a->nPoints, sizeof(struct kd_point),
                    _compPoints, axis, a->max_threads))
        return NULL;
    pivot = a->nPoints / 2;
    if ((node = kd_allocNode(a->points, pivot, a->min, a->max, axis,
                             a->dim)) == NULL)
        return NULL;

    arg[0] = _childArg(a, a->points, pivot, axis, node->location[axis], 0);
    arg[1] = _childArg(a, a->points + pivot, a->nPoints - pivot, axis,
                       node->location[axis], 1);
    if (!arg[0] || !arg[1]) {
        _freeChildArg(arg[0]);
        _freeChildArg(arg[1]);
        kd_freeNode(node, NULL);
        return NULL;
    }

    /*
     * Build each child in a new thread if allowed. A child whose
     * thread cannot be started is built in this thread.
     */
    for(i = 0; i < 2; i++)
        threaded[i] = a->max_threads > 1 &&
            !pthread_create(&threads[i], NULL, kd_doBuildTree, arg[i]);
    for(i = 0; i < 2; i++) {
        if (threaded[i]) {
            pthread_join(threads[i], &ret);
            child[i] = ret;
        } else {
            child[i] = kd_doBuildTree(arg[i]);
        }
        _freeChildArg(arg[i]);
    }

    node->left = child[0];
    node->right = child[1];
    if (!node->left || !node->right) {
        kd_destroyTree(node, a->destr);
        return NULL;
    }
    return node;
}




void
kd_freeNode(kdNode *node, void (*destr) (void *))
{
    if (!node)
	return;
    if (node->location)
        free(node->location);
    if (node->min)
        free(node->min);
    if (node->max)
        free(node->max);
    if (destr)
        destr(node->data);
    free(node);
    return;
}


struct kd_thread_data *
kd_buildArg(struct kd_point *points,
            unsigned long nPoints,
            void *(*constr) (void *),
            void (*destr) (void *), float *min,
            float *max, int depth, int max_threads, int dim)
{
    struct kd_thread_data *d;

    if ((d = kd_malloc(sizeof(kd_thread_data), "kd_thread_data")) == NULL)
        return NULL;
    d->points = points;
    d->nPoints = nPoints;
    d->constr = constr;
    d->destr = destr;
    d->min = min;
    d->max = max;
    d->depth = depth;
    d->max_threads = max_threads;
    d->dim = dim;
    return d;
}


struct kdNode *
kd_allocNode(struct kd_point *points, unsigned long pivot,
             float *min, float *max, int axis, int dim)
{
    struct kdNode *node;

    if ((node = kd_malloc(sizeof(kdNode), "kd_allocNode (node): ")) == NULL)
        return NULL;
    node->location = node->min = node->max = NULL;
    node->left = node->right = NULL;
    node->data = NULL;
    node->split = axis;
    if ((node->location = kd_malloc(sizeof(float) * dim,
                                    "kd_allocNode (node->location): ")) ==
        NULL) {
        kd_freeNode(node, NULL);
        return NULL;
    }
    memcpy(node->location, points[pivot].point, dim * sizeof(float));
    if ((node->min = kd_malloc(sizeof(float) * dim,
                               "kd_allocNode (node->min): ")) == NULL) {
        kd_freeNode(node, NULL);
        return NULL;
    }
    memcpy(node->min, min, dim * sizeof(float));
    if ((node->max = kd_malloc(sizeof(float) * dim,
                               "kd_allocNode (node->max): ")) == NULL) {
        kd_freeNode(node, NULL);
        return NULL;
    }
    memcpy(node->max, max, dim * sizeof(float));
    return node;
}

/*!
 * \brief free the kd-tree data structure, 
 *
 * \param node the root node of the tree to be destroyed
 *
 * \param *destr a pointer to the destructor function for the data container.
 *
 * \return This function does not return a value
 */
void
kd_destroyTree(struct kdNode *node, void (*destr) (void *))
{
    if (node == NULL)
        return;

    kd_destroyTree(node->left, destr);
    kd_destroyTree(node->right, destr);
    kd_freeNode(node, destr);
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
        if (!pqinsert(res, point)) {
            free(point);
            return 0;
        }
    }
    if (!kd_insertResTree(node->right, res))
        return 0;
    return 1;
}

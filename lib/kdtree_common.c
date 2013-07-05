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

void *
kd_doBuildTree(void *threadarg)
{

    int sortaxis;
    unsigned long pivot;
    float *tmpMinLeft, *tmpMaxLeft, *tmpMinRight, *tmpMaxRight;
    struct kdNode *node;
    struct kd_point *points;
    unsigned long nPoints;
    void *(*constr) (void *);
    void (*destr) (void *);
    float *min, *max;
    int depth, dim, max_threads;
    pthread_t threads[2];
    pthread_attr_t attr;
    struct kd_thread_data *argleft, *argright;
    struct kd_thread_data *my_data = (struct kd_thread_data *) threadarg;

    points = my_data->points;
    nPoints = my_data->nPoints;
    constr = my_data->constr;
    destr = my_data->destr;
    min = my_data->min;
    max = my_data->max;
    depth = my_data->depth;
    max_threads = my_data->max_threads;
    dim = my_data->dim;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    sortaxis = depth % dim;
    if (nPoints == 1) {
        if ((node = kd_allocNode(points, 0, min, max, sortaxis, dim)) == NULL)
            return NULL;
        if (constr)
            node->data = constr(points[0].data);
        return node;
    }

    /*
     * If this iteration is allowed to start more threads, we first
     * use them to parallelize the sorting 
     */
    pmergesort(points, nPoints, sizeof(struct kd_point), _compPoints,
               sortaxis, max_threads);
    pivot = nPoints / 2;
    if ((node = kd_allocNode(points, pivot, min, max, sortaxis, dim)) == NULL)
        return NULL;

    if ((tmpMinLeft = kd_malloc(sizeof(float) * dim,
                                "kd_buildTree (tmpMin): ")) == NULL) {
        kd_destroyTree(node, destr);
        return NULL;
    }
    memcpy(tmpMinLeft, min, dim * sizeof(float));
    if ((tmpMaxLeft = kd_malloc(sizeof(float) * dim,
                                "kd_buildTree (tmpMax): ")) == NULL) {
        kd_destroyTree(node, destr);
        free(tmpMinLeft);
        return NULL;
    }
    memcpy(tmpMaxLeft, max, dim * sizeof(float));
    tmpMaxLeft[sortaxis] = node->location[sortaxis];
    argleft =
        kd_buildArg(points, pivot, constr, destr, tmpMinLeft, tmpMaxLeft,
                    depth + 1, max_threads / 2, dim);
    if (!argleft) {
        kd_destroyTree(node, destr);
        free(tmpMinLeft);
        free(tmpMaxLeft);
        return NULL;
    }
    if (max_threads > 1) {
        pthread_create(&threads[0], &attr, kd_doBuildTree, (void *) argleft);
    } else {
        node->left = kd_doBuildTree((void *) argleft);
        free(argleft);
        free(tmpMinLeft);
        free(tmpMaxLeft);
        if (!node->left) {
            kd_destroyTree(node, destr);
            return NULL;
        }
    }

    if ((tmpMinRight = kd_malloc(sizeof(float) * dim,
                                 "kd_buildTree (tmpMin): ")) == NULL) {
        kd_destroyTree(node, destr);
        return NULL;
    }
    memcpy(tmpMinRight, min, dim * sizeof(float));
    if ((tmpMaxRight = kd_malloc(sizeof(float) * dim,
                                 "kd_buildTree (tmpMax): ")) == NULL) {
        kd_destroyTree(node, destr);
        free(tmpMinRight);
        return NULL;
    }

    memcpy(tmpMaxRight, max, dim * sizeof(float));
    tmpMinRight[sortaxis] = node->location[sortaxis];
    argright = kd_buildArg(&points[pivot], nPoints - pivot, constr, destr,
                           tmpMinRight, tmpMaxRight, depth + 1,
                           max_threads / 2, dim);
    if (!argright) {
        kd_destroyTree(node, destr);
        free(tmpMinRight);
        free(tmpMaxRight);
        return NULL;
    }
    if (max_threads > 1) {
        pthread_create(&threads[1], &attr, kd_doBuildTree, (void *) argright);
    } else {
        node->right = kd_doBuildTree((void *) argright);
        free(argright);
        free(tmpMinRight);
        free(tmpMaxRight);
        if (!node->right) {
            kd_destroyTree(node, destr);
            return NULL;
        }
    }
    if (max_threads > 1) {
        pthread_join(threads[0], (void *) (&node->left));
        free(argleft);
        free(tmpMinLeft);
        free(tmpMaxLeft);
        pthread_join(threads[1], (void *) (&node->right));
        free(argright);
        free(tmpMinRight);
        free(tmpMaxRight);
        if (!node->left || !node->right) {
            kd_destroyTree(node, destr);
            return NULL;
        }
    }
    return (void *) node;
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
    node->left = node->right = NULL;
    node->data = NULL;
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
        pqinsert(res, point);
    }
    if (!kd_insertResTree(node->right, res))
        return 0;
    return 1;
}

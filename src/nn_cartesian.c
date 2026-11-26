/* Nearest-neighbor searches in 3-d Euclidean space */

#include <stdio.h>
#include <stdlib.h>

#include <kdtree.h>

#define NPOINTS 1000
#define DIM 3

/* Copies the point's id into the tree */
static void *
copy_id(void *data)
{
    int *id = malloc(sizeof *id);

    if (id)
        *id = *(int *) data;
    return id;
}

int
main(void)
{
    static float coords[NPOINTS][DIM];
    static int ids[NPOINTS];
    struct kd_point points[NPOINTS];
    float min[DIM] = { 0, 0, 0 }, max[DIM] = { 1, 1, 1 };
    float query[DIM] = { 0.5, 0.5, 0.5 };
    float dist_sq;
    struct kdNode *tree, *nearest;
    struct pqueue *res;
    struct resItem *item;
    int i, j;

    /* Random points in the unit cube, each tagged with an id */
    for(i = 0; i < NPOINTS; i++) {
        for(j = 0; j < DIM; j++)
            coords[i][j] = rand() / (float) RAND_MAX;
        ids[i] = i;
        points[i].point = coords[i];
        points[i].data = &ids[i];
    }

    /* min and max bound all points. One thread builds the tree. */
    tree = kd_buildTree(points, NPOINTS, copy_id, free, min, max, DIM, 1);
    if (!tree)
        return EXIT_FAILURE;

    /* Search radius (squared) larger than the cube's diagonal, so a
       neighbor is always found. On return dist_sq holds the squared
       distance to the nearest neighbor. */
    dist_sq = DIM;
    nearest = kd_nearest(tree, query, &dist_sq, DIM);
    printf("Nearest neighbor: point %d at distance %.4f\n",
           *(int *) nearest->data, sqrt(dist_sq));

    /* The 5 nearest neighbors come back in a priority queue */
    dist_sq = DIM;
    if (!(res = kd_qnearest(tree, query, &dist_sq, 5, DIM)))
        return EXIT_FAILURE;

    /* pqremove_min hands out the items by increasing distance */
    while (pqremove_min(res, &item)) {
        printf("  point %d at distance %.4f\n",
               *(int *) item->node->data, sqrt(item->dist_sq));
        free(item);
    }
    free(res->d);
    free(res);

    kd_destroyTree(tree, free);
    return EXIT_SUCCESS;
}

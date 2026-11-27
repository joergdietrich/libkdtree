#include <stdio.h>
#include <stdlib.h>

#include <kdtree_internal.h>

#define NITEMS 100

int
main()
{
    struct kd_point point;
    struct pqueue *q;
    struct resItem *items[NITEMS], *item;
    float coords[2] = { 0.5, 0.5 };
    float min[2] = { 0, 0 }, max[2] = { 1, 1 };
    float last;
    int i;

    /* A tree needs at least one point */
    point.point = coords;
    point.data = NULL;
    if (kd_buildTree(&point, 0, NULL, NULL, min, max, 2, 1) != NULL ||
        kd_sph_buildTree(&point, 0, NULL, NULL, min, max, 1) != NULL) {
        fprintf(stderr, "Tree built from zero points\n");
        exit(EXIT_FAILURE);
    }

    /* A tree needs a bounding box */
    if (kd_buildTree(&point, 1, NULL, NULL, NULL, max, 2, 1) != NULL ||
        kd_sph_buildTree(&point, 1, NULL, NULL, min, NULL, 1) != NULL) {
        fprintf(stderr, "Tree built without a bounding box\n");
        exit(EXIT_FAILURE);
    }

    /* A queue initialized for 0 items still grows */
    if ((q = pqinit(NULL, 0)) == NULL)
        exit(EXIT_FAILURE);
    for(i = 0; i < NITEMS; i++) {
        items[i] = malloc(sizeof(struct resItem));
        items[i]->node = NULL;
        items[i]->dist_sq = (i * 37) % NITEMS;
        if (!pqinsert(q, items[i])) {
            fprintf(stderr, "Insert %d failed\n", i);
            exit(EXIT_FAILURE);
        }
    }
    last = -1;
    for(i = 0; i < NITEMS / 2; i++) {
        if (!pqremove_min(q, &item) || item->dist_sq < last) {
            fprintf(stderr, "Queue out of order at item %d\n", i);
            exit(EXIT_FAILURE);
        }
        last = item->dist_sq;
        free(item);
    }

    /* pqfree frees the items still in the queue */
    pqfree(q);
    pqfree(NULL);
    return 0;
}

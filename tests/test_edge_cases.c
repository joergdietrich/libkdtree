#include <stdio.h>
#include <stdlib.h>

#include <kdtree_internal.h>

#define NITEMS 100

/* Returns a queue initialized for 0 items, which must grow to hold
   the NITEMS items inserted in scrambled order */
static struct pqueue *
fill_queue(void)
{
    struct pqueue *q;
    struct resItem *item;
    int i;

    if ((q = pqinit(NULL, 0)) == NULL)
        exit(EXIT_FAILURE);
    for(i = 0; i < NITEMS; i++) {
        item = malloc(sizeof(struct resItem));
        item->node = NULL;
        item->dist_sq = (i * 37) % NITEMS;
        if (!pqinsert(q, item)) {
            fprintf(stderr, "Insert %d failed\n", i);
            exit(EXIT_FAILURE);
        }
    }
    return q;
}

int
main()
{
    struct kd_point point;
    struct pqueue *q;
    struct resItem *item;
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

    /* A queue initialized for 0 items still grows, and drains in order */
    q = fill_queue();
    last = -1;
    for(i = 0; i < NITEMS; i++) {
        if (!pqremove_min(q, &item) || item->dist_sq < last) {
            fprintf(stderr, "Queue out of order at item %d\n", i);
            exit(EXIT_FAILURE);
        }
        last = item->dist_sq;
        free(item);
    }
    if (pqremove_min(q, &item)) {
        fprintf(stderr, "Drained queue not empty\n");
        exit(EXIT_FAILURE);
    }
    pqfree(q);

    /* pqfree frees the items still in a partly drained queue */
    q = fill_queue();
    for(i = 0; i < NITEMS / 2; i++) {
        pqremove_max(q, &item);
        free(item);
    }
    pqfree(q);
    pqfree(NULL);
    return 0;
}

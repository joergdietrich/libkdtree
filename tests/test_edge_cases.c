#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <kdtree_internal.h>

#define NITEMS 100
#define NPOINTS 200

/* Returns a queue initialized for 0 items, which must grow to hold
   the NITEMS items inserted in scrambled order. Each growth must at
   least double the queue. */
static struct pqueue *
fill_queue(void)
{
    struct pqueue *q;
    struct resItem *item;
    uint32_t avail;
    int i;

    if ((q = kd_pqinit(NULL, 0)) == NULL)
        exit(EXIT_FAILURE);
    for(i = 0; i < NITEMS; i++) {
        if ((item = malloc(sizeof(struct resItem))) == NULL)
            exit(EXIT_FAILURE);
        item->node = NULL;
        item->dist_sq = (i * 37) % NITEMS;
        avail = q->avail;
        if (!kd_pqinsert(q, item)) {
            fprintf(stderr, "Insert %d failed\n", i);
            exit(EXIT_FAILURE);
        }
        if (q->avail != avail && q->avail < 2 * avail) {
            fprintf(stderr, "Queue grew from %u to only %u slots\n",
                    avail, q->avail);
            exit(EXIT_FAILURE);
        }
    }
    return q;
}

/* Checks that asking a tree of NPOINTS points for UINT_MAX nearest
   neighbors returns all of them */
static void
check_all_nearest(int spherical)
{
    static struct kd_point points[NPOINTS];
    static float coords[NPOINTS][2];
    float min[2] = { 0, -M_PI / 2 }, max[2] = { 2 * M_PI, M_PI / 2 };
    float p[2] = { 1, 0 }, range = 100;
    struct kdNode *tree;
    struct pqueue *res;
    int i;

    srand48(7);
    for(i = 0; i < NPOINTS; i++) {
        coords[i][0] = min[0] + drand48() * (max[0] - min[0]);
        coords[i][1] = min[1] + drand48() * (max[1] - min[1]);
        points[i].point = coords[i];
        points[i].data = NULL;
    }
    if (spherical) {
        tree = kd_sph_buildTree(points, NPOINTS, NULL, NULL, min, max, 1);
        res = tree ? kd_sph_qnearest(tree, p, &range, UINT_MAX) : NULL;
    } else {
        tree = kd_buildTree(points, NPOINTS, NULL, NULL, min, max, 2, 1);
        res = tree ? kd_qnearest(tree, p, &range, UINT_MAX, 2) : NULL;
    }
    if (!res || res->size - 1 != NPOINTS) {
        fprintf(stderr, "%s q-nearest search for UINT_MAX points found %d\n",
                spherical ? "Spherical" : "Cartesian",
                res ? (int) res->size - 1 : -1);
        exit(EXIT_FAILURE);
    }
    pqfree(res);
    kd_destroyTree(tree, NULL);
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

    /* q + 1 must not wrap around */
    check_all_nearest(0);
    check_all_nearest(1);

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

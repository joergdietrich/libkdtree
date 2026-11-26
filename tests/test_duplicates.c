#include <stdio.h>
#include <stdlib.h>

#include <kdtree_internal.h>

#define NPOINTS 20000
#define NQUERIES 200

/* Compares nearest-neighbor and range searches with brute force on
   points that share many coordinates. Returns 0 on success. */
static int
check(float (*coords)[2], const char *name, int nthreads)
{
    static struct kd_point points[NPOINTS];
    float min[2] = { 0, 0 }, max[2] = { 1, 1 };
    struct kdNode *tree, *nearest;
    struct pqueue *res;
    int i, k, inrange, failed = 0;
    uint32_t j;

    for(i = 0; i < NPOINTS; i++) {
        points[i].point = coords[i];
        points[i].data = NULL;
    }
    if ((tree = kd_buildTree(points, NPOINTS, NULL, NULL, min, max, 2,
                             nthreads)) == NULL) {
        fprintf(stderr, "%s: building the tree failed\n", name);
        return 1;
    }
    for(k = 0; k < NQUERIES && !failed; k++) {
        float q[2] = { drand48(), drand48() };
        float best = 4, dist_sq = 4, radius_sq = 0.01;

        inrange = 0;
        for(i = 0; i < NPOINTS; i++) {
            float d = kd_dist_sq(coords[i], q, 2);

            if (d < best)
                best = d;
            if (d < radius_sq)
                inrange++;
        }
        nearest = kd_nearest(tree, q, &dist_sq, 2);
        if (!nearest || kd_dist_sq(nearest->location, q, 2) != best) {
            fprintf(stderr, "%s, %d threads: wrong nearest neighbor\n",
                    name, nthreads);
            failed = 1;
            break;
        }
        if ((res = kd_range(tree, q, &radius_sq, 2, KD_UNORDERED)) == NULL) {
            failed = 1;
            break;
        }
        if ((int) res->size - 1 != inrange) {
            fprintf(stderr, "%s, %d threads: range search found %u of %d\n",
                    name, nthreads, res->size - 1, inrange);
            failed = 1;
        }
        for(j = 1; j < res->size; j++)
            free(res->d[j]);
        free(res->d);
        free(res);
    }
    kd_destroyTree(tree, NULL);
    return failed;
}

int
main()
{
    static float grid[NPOINTS][2], same[NPOINTS][2];
    int i, nthreads;

    srand48(5);
    for(i = 0; i < NPOINTS; i++) {
        /* 10 x 10 grid, about 200 points per grid node */
        grid[i][0] = (int) (drand48() * 10) / 10.;
        grid[i][1] = (int) (drand48() * 10) / 10.;
        same[i][0] = same[i][1] = 0.5;
    }
    for(nthreads = 1; nthreads <= 4; nthreads *= 4) {
        if (check(grid, "grid", nthreads) || check(same, "identical",
                                                   nthreads))
            exit(EXIT_FAILURE);
    }
    return 0;
}

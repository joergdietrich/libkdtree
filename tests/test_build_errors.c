#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include <kdtree.h>

#define NPOINTS 1000

/* The build calls copy() from several threads */
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static int ncopies, nfreed, fail_at;

/* Copies an int, failing on call number fail_at */
static void *
copy(void *data)
{
    int *p, n;

    pthread_mutex_lock(&lock);
    n = ++ncopies;
    pthread_mutex_unlock(&lock);
    if (n == fail_at || (p = malloc(sizeof(int))) == NULL)
        return NULL;
    *p = *(int *) data;
    return p;
}

static void
destroy(void *data)
{
    nfreed++;
    free(data);
}

int
main()
{
    static struct kd_point points[NPOINTS];
    static float coords[NPOINTS][2];
    static int ids[NPOINTS];
    float min[2] = { 0, 0 }, max[2] = { 1, 1 };
    struct kdNode *tree;
    int i, nthreads;

    srand48(3);
    for(i = 0; i < NPOINTS; i++) {
        coords[i][0] = drand48();
        coords[i][1] = drand48();
        ids[i] = i;
        points[i].point = coords[i];
        points[i].data = &ids[i];
    }

    /* A NaN coordinate makes the build fail */
    coords[NPOINTS / 2][1] = NAN;
    if (kd_buildTree(points, NPOINTS, NULL, NULL, min, max, 2, 1) != NULL ||
        kd_sph_buildTree(points, NPOINTS, NULL, NULL, min, max, 1) != NULL) {
        fprintf(stderr, "Tree built from a NaN coordinate\n");
        exit(EXIT_FAILURE);
    }
    coords[NPOINTS / 2][1] = 0.5;

    for(nthreads = 1; nthreads <= 4; nthreads *= 4) {
        /* A failing constr makes the build fail, and the copies made
           before it are freed */
        ncopies = nfreed = 0;
        fail_at = NPOINTS / 2;
        if (kd_buildTree(points, NPOINTS, copy, destroy, min, max, 2,
                         nthreads) != NULL) {
            fprintf(stderr, "Tree built although constr failed\n");
            exit(EXIT_FAILURE);
        }
        if (nfreed != ncopies - 1) {
            fprintf(stderr, "%d of %d data copies freed after a failed "
                    "build\n", nfreed, ncopies - 1);
            exit(EXIT_FAILURE);
        }
        if (nthreads == 1 && ncopies != fail_at) {
            fprintf(stderr, "constr called %d times after it failed\n",
                    ncopies - fail_at);
            exit(EXIT_FAILURE);
        }

        /* With a working constr every point gets its data */
        ncopies = nfreed = 0;
        fail_at = 0;
        if ((tree = kd_buildTree(points, NPOINTS, copy, destroy, min, max,
                                 2, nthreads)) == NULL) {
            fprintf(stderr, "Building the tree failed\n");
            exit(EXIT_FAILURE);
        }
        kd_destroyTree(tree, destroy);
        if (nfreed != NPOINTS) {
            fprintf(stderr, "%d of %d data copies freed\n", nfreed, NPOINTS);
            exit(EXIT_FAILURE);
        }
    }
    return 0;
}

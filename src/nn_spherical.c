/* Nearest-neighbor and range searches on the sphere */

#include <stdio.h>
#include <stdlib.h>

#include <kdtree.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NPOINTS 1000
#define DEG (M_PI / 180)

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
    static float coords[NPOINTS][2];
    static int ids[NPOINTS];
    struct kd_point points[NPOINTS];

    /* Longitude (RA) in [0, 2 pi], latitude (Dec) in [-pi/2, pi/2] */
    float min[2] = { 0, -M_PI / 2 }, max[2] = { 2 * M_PI, M_PI / 2 };
    float query[2] = { 120 * DEG, -30 * DEG };
    float dist_sq;
    struct kdNode *tree, *nearest;
    struct pqueue *res;
    unsigned int i;

    /* Random points spread evenly over the sphere */
    for(i = 0; i < NPOINTS; i++) {
        coords[i][0] = 2 * M_PI * rand() / RAND_MAX;
        coords[i][1] = asin(2.0 * rand() / RAND_MAX - 1);
        ids[i] = i;
        points[i].point = coords[i];
        points[i].data = &ids[i];
    }

    tree = kd_sph_buildTree(points, NPOINTS, copy_id, free, min, max, 1);
    if (!tree)
        return EXIT_FAILURE;

    /* Distances are great-circle angles in radians. No two points are
       more than pi apart, so this radius always finds a neighbor. */
    dist_sq = M_PI * M_PI;
    nearest = kd_sph_nearest(tree, query, &dist_sq);
    printf("Nearest neighbor: point %d, %.3f deg away\n",
           *(int *) nearest->data, sqrt(dist_sq) / DEG);

    /* All points within 10 degrees, in no particular order */
    dist_sq = (10 * DEG) * (10 * DEG);
    if (!(res = kd_sph_range(tree, query, &dist_sq, KD_UNORDERED)))
        return EXIT_FAILURE;

    /* Items sit at res->d[1] to res->d[res->size - 1] */
    printf("%u points within 10 deg\n", res->size - 1);
    for(i = 1; i < res->size; i++)
        free(res->d[i]);
    free(res->d);
    free(res);

    kd_destroyTree(tree, free);
    return EXIT_SUCCESS;
}

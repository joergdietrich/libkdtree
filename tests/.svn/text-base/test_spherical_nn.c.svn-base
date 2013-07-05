/*! \file test_spherical_nn.c
 * \brief compare kdTree NN search with naive NN.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <kdtree.h>

#include "kdtree_test.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

void
naive_sph_nn(struct kd_point *pointlist, unsigned long npoints,
             float *point, float *nn, float *dist_sq)
{
    unsigned long i;
    float dsq;

    *dist_sq = kd_sph_dist_sq(pointlist[0].point, point);
    memcpy(nn, pointlist[0].point, 2 * sizeof(float));

    for(i = 1; i < npoints; i++) {
        dsq = kd_sph_dist_sq(pointlist[i].point, point);
        if (dsq < *dist_sq) {
            *dist_sq = dsq;
            memcpy(nn, pointlist[i].point, 2 * sizeof(float));
        }
    }
}


int
main(int argc, char **argv)
{
    unsigned int npoints;
    unsigned int nthreads;
    struct kd_point *pointlist;
    struct kdNode *kdTree;
    struct kdNode *nearest = NULL;
    float min[2], max[2];
    float point[2];
    float nn[2];
    float range, nnrange;

    if (argc != 3) {
        fprintf(stderr,
                "Usage: %s <no of points> <no of threads>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    npoints = atoi(argv[1]);
    nthreads = atoi(argv[2]);
    /*
     * Cover the entire sphere with points
     */
    min[0] = 0;
    min[1] = -M_PI / 2.;
    max[0] = 2 * M_PI;
    max[1] = M_PI / 2.;
    if (!(pointlist = create_pointlist(npoints, min, max, 2)))
        exit(EXIT_FAILURE);

    /*
     * Now we build the kd-tree from the pointlist 
     */
    if ((kdTree = kd_sph_buildTree(pointlist, npoints, NULL, NULL, min, max,
                               nthreads)) == NULL) {
        fprintf(stderr, "Error building kd-tree\n");
        exit(EXIT_FAILURE);
    }

    /*
     * Search nearest neighbors on a grid 
     */
    for(point[0] = 0; point[0] < 2 * M_PI; point[0] += 0.5) {
        for(point[1] = -M_PI / 2; point[1] < M_PI / 2; point[1] += 0.5) {
            range = nnrange = 4 * M_PI * M_PI;
            nearest = kd_sph_nearest(kdTree, point, &range);
            naive_sph_nn(pointlist, npoints, point, nn, &nnrange);
            if (!points_eq(nearest, nn, 2) || fabs(range - nnrange) > EPS) {
		fprintf(stderr, "Nearest neighbor to point (%.5f, %.5f):\n",
			point[0], point[1]);
		fprintf(stderr, "kd-tree: (%.5f, % .5f), dist_sq = %.5f\n",
			nearest->location[0], nearest->location[1], range);
		fprintf(stderr, "naive:   (%.5f, % .5f), dist_sq = %.5f\n",
			nn[0], nn[1], nnrange);
                exit(EXIT_FAILURE);
	    }
        }
    }
    return 0;
}

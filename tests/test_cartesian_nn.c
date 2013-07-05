/*! \file test_cartesian_nn.c
 * \brief compare kdTree NN search with naive NN.
 */
#include <stdlib.h>
#include <string.h>

#include <kdtree.h>

#include "kdtree_test.h"

#define DIM 2

void
naive_nn(struct kd_point *pointlist, unsigned long npoints, 
	 unsigned short dim, float *point, float *nn, float *dist_sq)
{
    unsigned long i;
    float dsq;
    
    *dist_sq = kd_dist_sq(pointlist[0].point, point, dim);
    memcpy(nn, pointlist[0].point, dim*sizeof(float));
    
    for(i=1; i<npoints; i++) {
	dsq = kd_dist_sq(pointlist[i].point, point, dim);
	if (dsq < *dist_sq) {
	    *dist_sq = dsq;
	    memcpy(nn, pointlist[i].point, dim*sizeof(float));
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
    float min[DIM], max[DIM];
    float point[DIM];
    float nn[DIM];
    float range, nnrange;
    unsigned short dim = DIM;

    if (argc != 3) {
	fprintf(stderr,
                "Usage: %s <no of points> <no of threads>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    npoints = atoi(argv[1]);
    nthreads = atoi(argv[2]);
    /*
     * All point will lie within the unit square 
     */
    min[0] = min[1] = 0;
    max[0] = max[1] = 1;
    if (!(pointlist = create_pointlist(npoints, min, max, dim)))
	exit(EXIT_FAILURE);
    
    /*
     * Now we build the kd-tree from the pointlist 
     */
    if ((kdTree = kd_buildTree(pointlist, npoints, NULL, NULL, min, max, dim, 
			       nthreads)) == NULL) {
        fprintf(stderr, "Error building kd-tree\n");
        exit(EXIT_FAILURE);
    }

    /* Search nearest neighbors on a 10 x 10 grid */
    for (point[0]=0; point[0]<1; point[0]+=0.1) {
	for(point[1]=0; point[1]<1; point[1]+=0.1) {
	    range = nnrange = 2;
	    nearest = kd_nearest(kdTree, point, &range, dim);
	    naive_nn(pointlist, npoints, dim, point, nn, &nnrange);
	    if (!points_eq(nearest, nn, dim) || fabs(range - nnrange) > EPS) {
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

#include <stdlib.h>
#include <string.h>

#include <kdtree.h>
#include <pqueue.h>

#include "kdtree_test.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif
#define NTEST 5

struct pqueue *
naive_sph_range(struct kd_point *pointlist, unsigned int npoints, float *p, 
		float *range)
{
    struct pqueue *res;
    struct resItem *point;
    float dsq;
    unsigned int i;

    if ((res = pqinit(NULL, 1)) == NULL) {
	fprintf(stderr, "Could not allocate memory for result queue\n.");
	return NULL;
    }
    for(i=0; i<npoints; i++) {
	if ((dsq = kd_sph_dist_sq(pointlist[i].point, p)) < *range) {
	    if ((point = kd_malloc(sizeof(struct resItem), "naive_qnearest: "))
		== NULL)
		return NULL;
	    if ((point->node = kd_allocNode(pointlist, i, 
					    /* 2 dummies */
					    pointlist[i].point, 
					    pointlist[i].point,
					    -1, 2)) == NULL)
		return NULL;
	    point->dist_sq = dsq;
	    pqinsert(res, point);
	}
    }
    return res;
}

int
main(int argc, char **argv)
{
    unsigned int npoints;
    unsigned int nthreads;
    struct kd_point *pointlist;
    struct kdNode *kdTree;
    float min[2], max[2];
    float point[2];
    struct pqueue *result, *nresult;
    unsigned int i;
    float ranges[NTEST] = {0.05, 0.1, 0.25, 0.5, 1};

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
    if ((kdTree = kd_buildTree(pointlist, npoints, NULL, NULL, min, max, 2, 
			       nthreads)) == NULL) {
        fprintf(stderr, "Error building kd-tree\n");
        exit(EXIT_FAILURE);
    }

    /* Search q-nearest neighbors on a grid */
    for (point[0]=min[0]; point[0]<max[0]; point[0]+=0.5) {
	for(point[1]=min[1]; point[1]<max[1]; point[1]+=0.5) {
	    for(i=0; i<NTEST; i++) {
		if ((result = kd_sph_range(kdTree, point, &ranges[i],
					   KD_ORDERED)) == NULL) {
		    fprintf(stderr, 
			    "kd_sph_range returned NULL for range search around (%.5f, %.5f), radius sq. %.5f with %d threads.\n",
			    point[0], point[1], ranges[i], nthreads);
		    exit(EXIT_FAILURE);
		}
		if ((nresult = naive_sph_range(pointlist, npoints, point, 
					       &ranges[i])) == NULL) {
		    fprintf(stderr, 
			    "naive_sph_range returned NULL for range search around (%.5f, %.5f), radius sq. %.5f with %d threads.\n",
			    point[0], point[1], ranges[i], nthreads);
		    exit(EXIT_FAILURE);
		}
		if (!sorted_queues_eq(result, nresult, 2)) {
		    fprintf(stderr, "%d points, %d threads, range search radius sq. %.5f.\n",
			    npoints, nthreads, ranges[i]);
		    fprintf(stderr, "Searching around point (%.5f, %.5f).\n",
			    point[0], point[1]);
		    exit(EXIT_FAILURE);
		}
	    }
	}
    }
    
    return 0;
}


#include <stdlib.h>
#include <string.h>

#include <kdtree.h>
#include <pqueue.h>

#include "kdtree_test.h"

#define DIM 2
#define NTEST 5

struct pqueue *
naive_range(struct kd_point *pointlist, unsigned int npoints, 
	    unsigned short dim, float *p, float *range)
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
	if ((dsq = kd_dist_sq(pointlist[i].point, p, dim)) < *range) {
	    if ((point = kd_malloc(sizeof(struct resItem), "naive_qnearest: "))
		== NULL)
		return NULL;
	    if ((point->node = kd_allocNode(pointlist, i, 
					    /* 2 dummies */
					    pointlist[i].point, 
					    pointlist[i].point,
					    -1, dim)) == NULL)
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
    float min[DIM], max[DIM];
    float point[DIM];
    unsigned short dim = DIM;
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

    /* Range searches on a 10 x 10 grid */
    for (point[0]=0; point[0]<1; point[0]+=0.1) {
	for(point[1]=0; point[1]<1; point[1]+=0.1) {
	    for(i=0; i<NTEST; i++) {
		if ((result = kd_range(kdTree, point, &ranges[i], dim,
					  KD_ORDERED)) == NULL) {
		    fprintf(stderr, 
			    "kd_range returned NULL for range search around (%.5f, %.5f), radius sq. %.5f with %d threads.\n",
			    point[0], point[1], ranges[i], nthreads);
		    exit(EXIT_FAILURE);
		}
		if ((nresult = naive_range(pointlist, npoints, dim, point, 
					   &ranges[i])) == NULL) {
		    fprintf(stderr, 
			    "naive_range returned NULL for range search around (%.5f, %.5f), radius sq. %.5f with %d threads.\n",
			    point[0], point[1], ranges[i], nthreads);
		    exit(EXIT_FAILURE);
		}
		if (!sorted_queues_eq(result, nresult, DIM)) {
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


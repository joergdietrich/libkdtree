#include <stdlib.h>
#include <string.h>

#include <kdtree.h>
#include <pqueue.h>

#include "kdtree_test.h"

#define DIM 2
#define NTEST 5


struct pqueue *
naive_qnearest(struct kd_point *pointlist, unsigned int npoints, 
	       unsigned short dim, float *p, float *range, 
	       unsigned int q)
{
    struct pqueue *res;
    struct resItem *point, *item;
    float dsq;
    unsigned int i;

    if ((res = pqinit(NULL, q + 2)) == NULL) {
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
	if (res->size > q + 1) {
	    pqremove_max(res, &item);
	    free(item);
	    if (res->size > 1) {
		/*
		 * Only inspect the queue if there are items left 
		 */
		pqpeek_max(res, &item);
		*range = item->dist_sq;
	    } else {
		/* Nothing found */
		*range = 0;
	    }
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
    float range, nrange;
    unsigned short dim = DIM;
    struct pqueue *result, *nresult;
    unsigned int i;
    int qnearest[NTEST] = {5, 10, 20, 50, 100};

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

    /* Search q-nearest neighbors on a 10 x 10 grid */
    for (point[0]=0; point[0]<1; point[0]+=0.1) {
	for(point[1]=0; point[1]<1; point[1]+=0.1) {
	    for(i=0; i<NTEST; i++) {
		range = nrange = 4;
		if ((result = kd_qnearest(kdTree, point, &range, qnearest[i], 
					  dim)) == NULL) {
		    fprintf(stderr, 
			    "kd_qnearest returned NULL for %d-NN search around (%.5f, %.5f) with %d threads.\n",
			    qnearest[i], point[0], point[1], nthreads);
		    exit(EXIT_FAILURE);
		}
		if ((nresult = naive_qnearest(pointlist, npoints, dim, point, 
					      &nrange, qnearest[i])) == NULL) {
		    fprintf(stderr, 
			    "naive_qnearest returned NULL for %d-NN search around (%.5f, %.5f) with %d threads.\n",
			    qnearest[i], point[0], point[1], nthreads);
		    exit(EXIT_FAILURE);
		}
		if (!sorted_queues_eq(result, nresult, DIM)) {
		    fprintf(stderr, "%d points, %d threads, %d q-nearest neighbors.\n",
			    npoints, nthreads, qnearest[i]);
		    exit(EXIT_FAILURE);
		}
	    }
	}
    }
    return 0;
}

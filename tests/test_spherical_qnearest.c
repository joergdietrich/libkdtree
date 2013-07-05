#include <math.h>
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
naive_sph_qnearest(struct kd_point *pointlist, unsigned int npoints, 
		   float *p, float *range, unsigned int q)
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
    float min[2], max[2];
    float point[2];
    float range, nrange;
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
		range = nrange = 4 * M_PI * M_PI;;
		if ((result = kd_sph_qnearest(kdTree, point, &range, 
					      qnearest[i])) == NULL) {
		    fprintf(stderr, 
			    "kd_sph_qnearest returned NULL for %d-NN search around (%.5f, %.5f) with %d threads.\n",
			    qnearest[i], point[0], point[1], nthreads);
		    exit(EXIT_FAILURE);
		}
		if ((nresult = naive_sph_qnearest(pointlist, npoints, point, 
						  &nrange, qnearest[i])) 
		    == NULL) {
		    fprintf(stderr, 
			    "naive_sph_qnearest returned NULL for %d-NN search around (%.5f, %.5f) with %d threads.\n",
			    qnearest[i], point[0], point[1], nthreads);
		    exit(EXIT_FAILURE);
		}
		if (!sorted_queues_eq(result, nresult, 2)) {
		    fprintf(stderr, "%d points, %d threads, %d q-nearest neighbors.\n",
			    npoints, nthreads, qnearest[i]);
		    exit(EXIT_FAILURE);
		}
	    }
	}
    }
    return 0;
}

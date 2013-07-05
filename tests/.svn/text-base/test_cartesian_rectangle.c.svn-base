#include <stdlib.h>
#include <string.h>

#include <kdtree.h>
#include <pqueue.h>

#include "kdtree_test.h"

#define DIM 2
#define NTEST 5

int
main(int argc, char **argv)
{
    unsigned int npoints;
    unsigned int nthreads;
    struct kd_point *pointlist;
    struct kdNode *kdTree;
    float min[DIM], max[DIM];
    float point[DIM];
    float rect_min[DIM], rect_max[DIM];
    unsigned short dim = DIM;
    struct pqueue *result, *nresult;
    unsigned int i, j;
    float dx[NTEST] = {0.05, 0.1, 0.25, 0.5, 0.7};

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
    for (point[0]=0; point[0]<1; point[0]+=0.1) {
	for(point[1]=0; point[1]<1; point[1]+=0.1) {
	    for(i=0; i<NTEST; i++) {
		for(j=0; j<NTEST; j++) {
		    rect_min[0] = point[0] - dx[i];
		    rect_min[1] = point[1] - dx[j];
		    rect_max[0] = point[0] + dx[i];
		    rect_max[1] = point[1] + dx[j];
		    if ((result = kd_ortRangeSearch(kdTree, rect_min, 
						    rect_max, dim)) == NULL) {
			fprintf(stderr, 
				"Orthogonal range search returned NULL.\n");
			exit(EXIT_FAILURE);
		    }
		    if ((nresult = naive_ortRangeSearch(pointlist, npoints,
							rect_min, rect_max, 
							dim)) 
			== NULL) {
			fprintf(stderr, 
				"naive orthogonal range search returned NULL.\n");
			exit(EXIT_FAILURE);
		    }
		    if (!unsorted_queues_eq(result, nresult, DIM)) {
			fprintf(stderr, "%d points, %d threads, \n",
				npoints, nthreads);
			fprintf(stderr, 
				"rectangle (% .5f, % .5f) (% .5f, % .5f)\n",
				rect_min[0], rect_min[1], rect_max[0], 
				rect_max[1]);
			exit(EXIT_FAILURE);
		    }
		}
	    }
	}
    }
    return 0;
}

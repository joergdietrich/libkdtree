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


int
main(int argc, char **argv)
{
    unsigned int npoints;
    unsigned int nthreads;
    struct kd_point *pointlist;
    struct kdNode *kdTree;
    float min[2], max[2];
    float point[2];
    float rect_min[2], rect_max[2];
    struct pqueue *result, *nresult;
    unsigned int i, j;
    float dx[NTEST] = {0.1, M_PI/4, M_PI/2, M_PI, M_PI};
    float dy[NTEST] = {0.05, 0.1, 0.25, M_PI/4, M_PI/2};

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
    for (point[0]=0; point[0] < 2*M_PI; point[0]+=0.5) {
	for(point[1]=-M_PI/2; point[1]<M_PI/2; point[1]+=0.5) {
	    for(i=0; i<NTEST; i++) {
		for(j=0; j<NTEST; j++) {
		    rect_min[0] = point[0] - dx[i];
		    rect_min[1] = point[1] - dy[j];
		    rect_max[0] = point[0] + dx[i];
		    rect_max[1] = point[1] + dy[j];
		    /* Ensure orthogonal ranges don't go outside the boundary 
		     */
		    if (rect_min[0] < 0)
			rect_min[0] = 0;
		    if (rect_max[0] > 2 * M_PI)
			rect_max[0] = 2 * M_PI;
		    /* Ensure orthogonal ranges don't go over the poles */
		    if (rect_min[1] < -M_PI / 2)
			rect_min[1] =  -M_PI / 2;
		    if (rect_max[1] > M_PI / 2) 
			rect_max[1] = M_PI / 2.;
		    if ((result = kd_sph_ortRangeSearch(kdTree, rect_min, 
							rect_max)) == NULL) {
			fprintf(stderr, 
				"Orthogonal range search returned NULL.\n");
			exit(EXIT_FAILURE);
		    }
		    if ((nresult = naive_sph_ortRangeSearch(pointlist, npoints,
							    rect_min, 
							    rect_max, 2)) 
			== NULL) {
			fprintf(stderr, 
				"naive orthogonal range search returned NULL.\n");
			exit(EXIT_FAILURE);
		    }
		    if (!unsorted_queues_eq(result, nresult, 2)) {
			fprintf(stderr, "%d points, %d threads, \n",
				npoints, nthreads);
			fprintf(stderr, 
				"rectangle (% .5f, % .5f) (% .5f, % .5f)\n",
				rect_min[0], rect_min[1], rect_max[0], 
				rect_max[1]);
			fprintf(stderr, "i: %d\tj: %d\n", i, j);
			exit(EXIT_FAILURE);
		    }
		}
	    }
	}
    }
    return 0;
}

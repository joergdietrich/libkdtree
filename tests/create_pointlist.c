#include <stdlib.h>

#include <kdtree.h>

#include "kdtree_test.h"


struct kd_point *
create_pointlist(unsigned long npoints, float *min, float *max, 
		 unsigned short dim)
{
    unsigned long i;
    unsigned short j;
    struct kd_point *pointlist;
    
    srand48(11);
    pointlist = malloc(npoints * sizeof(struct kd_point));
    if (!pointlist)
	return NULL;

    for(i = 0; i < npoints; i++) {
        pointlist[i].point = malloc(dim * sizeof(float));
	if (!pointlist[i].point)
	    return NULL;
        pointlist[i].data = NULL;
        for(j = 0; j < dim; j++)
            pointlist[i].point[j] = drand48() * (max[j] - min[j]) + min[j];
    }
    return pointlist;
}

#include <kdtree_internal.h>

#include "kdtree_test.h"

int 
is_pointInRect(struct kd_point point, float *rect_min, float *rect_max,
	       unsigned short dim)
{
    unsigned short i;
    for(i=0; i<dim; i++) {
	if (point.point[i] < rect_min[i] || point.point[i] > rect_max[i])
	    return 0;
    }
    return 1;
}


struct pqueue *
naive_ortRangeSearch(struct kd_point *pointlist, unsigned int npoints,
		     float *rect_min, float *rect_max, unsigned short dim)
{
    struct pqueue *res;
    struct resItem *point;
    unsigned int i;

    if ((res = pqinit(NULL, 1)) == NULL ) {
	fprintf(stderr, "Could not allocate memory for result queue\n.");
	return NULL;
    }
    for(i=0; i<npoints; i++) {
	if (is_pointInRect(pointlist[i], rect_min, rect_max, dim)) {
	    if ((point = naive_item(pointlist[i].point, -1)) == NULL)
		return NULL;
            pqinsert(res, point);
	}
    }
    return res;
}


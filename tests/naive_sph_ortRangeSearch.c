#include <math.h>

#include <kdtree.h>
#include <pqueue.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

int 
sph_is_pointInRect(struct kd_point point, float *rect_min, float *rect_max,
		   unsigned short dim)
{
    if (point.point[1] < rect_min[1] || point.point[1] > rect_max[1])
	return 0;

    if (rect_min[0] < 0 && point.point[0] > M_PI) {
	if (point.point[0] - 2 * M_PI < rect_min[0] ||
	    point.point[0] - 2 * M_PI > rect_max[0])
	    return 0;
    } else {
	if (point.point[0] < rect_min[0] || point.point[0] > rect_max[0])
	    return 0;
    }

    return 1;
}


struct pqueue *
naive_sph_ortRangeSearch(struct kd_point *pointlist, unsigned int npoints,
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
	if (sph_is_pointInRect(pointlist[i], rect_min, rect_max, dim)) {
	    if ((point = kd_malloc(sizeof(struct resItem), 
				   "naive_ortRangeSearch: ")) == NULL)
                return NULL;
            if ((point->node = kd_allocNode(pointlist, i, 
                                            /* 2 dummies */
                                            pointlist[i].point, 
                                            pointlist[i].point,
                                            -1, dim)) == NULL)
                return NULL;
	    point->dist_sq = -1;
            pqinsert(res, point);
	}
    }
    return res;
}


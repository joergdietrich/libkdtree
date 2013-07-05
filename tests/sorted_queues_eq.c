#include <kdtree.h>

#include "kdtree_test.h"

int
sorted_queues_eq(struct pqueue *q1, struct pqueue *q2, unsigned short dim)
{
    struct resItem *p1, *p2;

    if (q1->size != q2->size) {
	fprintf(stderr, "Queues have different sizes: %d != %d\n", 
		q1->size, q2->size);
	return 0;
    }
    while(pqremove_min(q1, &p1)) {
	pqremove_min(q2, &p2);
	if (kd_dist_sq(p1->node->location, p2->node->location, dim) > EPS ||
	    fabs(p1->dist_sq - p2->dist_sq) > EPS) {
	    fprintf(stderr, "Points are not equal:\n");
	    fprintf(stderr, "Queue 1: (%.5f, %.5f) dist_sq: %.5f\n", 
		    p1->node->location[0], p1->node->location[1], p1->dist_sq);
	    fprintf(stderr, "Queue 2: (%.5f, %.5f) dist_sq: %.5f\n", 
		    p2->node->location[0], p2->node->location[1], p2->dist_sq);
	    return 0;
	}
	free(p1);
	free(p2);
    }
    return 1;
}

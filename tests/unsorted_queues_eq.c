#include <kdtree.h>

#include "kdtree_test.h"

int 
unsorted_queues_eq(struct pqueue *p1, struct pqueue *p2, unsigned short dim)
{
    uint32_t i;
    float *ref;
    struct pqueue *s1, *s2;

    /* Reference point at the origin from which we compute dummy
       distances to bring the queues into a sorted order. */
    ref = calloc(dim, sizeof(float));
    s1 = pqinit(NULL, 1);
    for(i=1; i<p1->size; i++) {
	p1->d[i]->dist_sq = kd_dist_sq(p1->d[i]->node->location, ref, dim);
	pqinsert(s1, p1->d[i]);
    }
    s2 = pqinit(NULL, 1);
    for(i=1; i<p2->size; i++) {
	p2->d[i]->dist_sq = kd_dist_sq(p2->d[i]->node->location, ref, dim);
	pqinsert(s2, p2->d[i]);
    }
    return sorted_queues_eq(s1, s2, dim);
}

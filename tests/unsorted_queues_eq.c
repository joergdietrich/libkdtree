#include <kdtree_internal.h>

#include "kdtree_test.h"

/* Inserts the items of p into s, with their distance from ref as the
   priority. Returns 0 if the queue could not grow. */
static int
resort(struct pqueue *p, struct pqueue *s, float *ref, unsigned short dim)
{
    uint32_t i;

    for(i = 1; i < p->size; i++) {
	p->d[i]->dist_sq = kd_dist_sq(p->d[i]->node->location, ref, dim);
	if (!kd_pqinsert(s, p->d[i]))
	    return 0;
    }
    return 1;
}

int
unsorted_queues_eq(struct pqueue *p1, struct pqueue *p2, unsigned short dim)
{
    float *ref;
    struct pqueue *s1, *s2;
    int ok;

    /* Reference point at the origin from which we compute dummy
       distances to bring the queues into a sorted order. */
    ref = calloc(dim, sizeof(float));
    s1 = kd_pqinit(NULL, 1);
    s2 = kd_pqinit(NULL, 1);
    ok = ref && s1 && s2 && resort(p1, s1, ref, dim) &&
	resort(p2, s2, ref, dim);
    free(ref);

    /* s1 and s2 only borrow the items. On failure p1 and p2 keep them,
       on success s1 and s2 take them over. */
    if (!ok) {
	fprintf(stderr, "Could not sort the queues for comparison\n");
	if (s1)
	    s1->size = 1;
	if (s2)
	    s2->size = 1;
	pqfree(s1);
	pqfree(s2);
	pqfree(p1);
	pqfree(p2);
	return 0;
    }
    p1->size = p2->size = 1;
    pqfree(p1);
    pqfree(p2);
    return sorted_queues_eq(s1, s2, dim);
}

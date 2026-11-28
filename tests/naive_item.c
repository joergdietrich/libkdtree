#include <stdlib.h>

#include <kdtree_internal.h>

#include "kdtree_test.h"

/* The item comes first, so free() on the item releases both. */
struct naive_item {
    struct resItem item;
    struct kdNode node;
};

struct resItem *
naive_item(float *location, float dist_sq)
{
    struct naive_item *p;

    if ((p = calloc(1, sizeof(struct naive_item))) == NULL)
	return NULL;
    p->node.location = location;
    p->item.node = &p->node;
    p->item.dist_sq = dist_sq;
    return &p->item;
}

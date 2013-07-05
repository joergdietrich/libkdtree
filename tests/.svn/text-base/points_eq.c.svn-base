#include <math.h>

#include <kdtree.h>

#include "kdtree_test.h"

/*! \brief Compare two points for equality within tolerance EPS.
 * \return 1 if equal, 0 if not.
 */
int
points_eq(struct kdNode *nkd, float *nn, unsigned short dim)
{
    unsigned short i;
    for(i=0; i<dim; i++) {
	if (fabs(nkd->location[i] - nn[i]) > EPS)
	    return 0;
    }
    return 1;
}


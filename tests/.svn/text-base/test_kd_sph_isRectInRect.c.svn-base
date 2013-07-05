#include <math.h>
#include <stdio.h>

#include <kdtree.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define NTESTS (3)

int
main()
{
    float minList[NTESTS][2] = { /* small rectangle, no meridian */
	{6./9. * M_PI, 5./18. * M_PI},
	{6./9. * M_PI, 5./18. * M_PI},
	/* small rectangle, meridian */
	{-2./36 * M_PI, -1./18. * M_PI}};
    float maxList[NTESTS][2] = {/* small rectangle, no meridian */
	{7./9. * M_PI, 7./18. * M_PI},
	{7./9. * M_PI, 7./18. * M_PI},
	/* small rectangle, meridian */
	{1./36. * M_PI, 2./18. * M_PI}};

    float minListTest[NTESTS][2] = {/* small rectangle contained */
	{6.2/9. * M_PI, 5.5/18. * M_PI},
	/* small rectangle, no meridian, uncontained (left) */
	{5.9/9. * M_PI, 5.5/18. * M_PI},
    	/* small rectangle, meridian, contained */
	{-1./36. * M_PI, 0},
    	/* small rectangle, meridian, uncontained (bottom) */
	{-1./36. * M_PI, -2./18. * M_PI}};

    float maxListTest[NTESTS][2] = {/* small rectangle contained */
	{6.8/9. * M_PI, 6.5/18. * M_PI},
	/* small rectangle, no meridian, uncontained (left) */
	{6.8/9. * M_PI, 6.5/18. * M_PI},
    	/* small rectangle, meridian, contained */
	{0.5/36. * M_PI, 1./18. * M_PI},
    	/* small rectangle, meridian, uncontained (bottom) */
	{0.5/36. * M_PI, 1./18. * M_PI}};

    short int truthList[NTESTS] = {
	/* small rectangle, no meridian, contained */
	1,
	/* small rectangle, no meridian, uncontained (left) */
	0,
    	/* small rectangle, meridian, contained */
	1,
    	/* small rectangle, meridian, uncontained (bottom) */
	0};


    int i;
    kd_point *point;
    kdNode *node;

    point = kd_malloc(sizeof(kd_point), "kd_point");
    point->point = kd_malloc(2 * sizeof(float), "point position");
    point->data = NULL;
    for(i = 0; i < NTESTS; i++) {
	point->point[0] = minListTest[i][0];
	point->point[1] = minListTest[i][1];
	node = kd_allocNode(point, 0, minListTest[i], maxListTest[i], 0, 2);
	if (truthList[i] != kd_sph_isRectInRect(node, minList[i],
						maxList[i])) {
	    fprintf(stderr, "Rectangle (%.3f, %.3f) (%.3f, %.3f)\n",
		    minList[i][0], minList[i][1], maxList[i][0], 
		    maxList[i][1]);
	    fprintf(stderr, "Test rectangle (%.3f, %.3f) (%.3f, %.3f)\n",
		    minListTest[i][0], minListTest[i][1], maxListTest[i][0], 
		    maxListTest[i][1]);
	    fprintf(stderr, "Expected truth value: %d\n", truthList[i]);
	    return 1;
	}
    }

    return 0;
}

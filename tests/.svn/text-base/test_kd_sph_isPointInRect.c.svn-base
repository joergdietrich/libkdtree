#include <math.h>
#include <stdio.h>

#include <kdtree.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define NTESTS (15)

int
main()
{
    float minList[NTESTS][2] = { /* small rectangle, no meridian */
	{2./3. * M_PI, 5./18. * M_PI},
	{2./3. * M_PI, 5./18. * M_PI},
	/* small rectangle, meridian */
	{-1./18 * M_PI, -1./18. * M_PI},
	{-1./18 * M_PI, -1./18. * M_PI},
	{-1./18 * M_PI, -1./18. * M_PI},
	/* long rectangle, no meridian */
	{1./18. * M_PI, -1./6. * M_PI},
	{1./18. * M_PI, -1./6. * M_PI},
	{1./18. * M_PI, -1./6. * M_PI},
	/* long rectangle, meridian */
	{-17./18. * M_PI, 1./9. * M_PI},
	{-17./18. * M_PI, 1./9. * M_PI},
	{-17./18. * M_PI, 1./9. * M_PI},
	/* rectangle around */
	{0, -1/18. * M_PI},
	{0, -1/18. * M_PI},
	/* triangle at north pole, over meridian */
	{-1./9. * M_PI, 6./18. * M_PI},
	{-1./9. * M_PI, 6./18. * M_PI}};
    float maxList[NTESTS][2] = {/* small rectangle, no meridian */
	{7./9. * M_PI, 7./18. * M_PI},
	{7./9. * M_PI, 7./18. * M_PI},
	/* small rectangle, meridian */
	{1./36. * M_PI, 1./9. * M_PI},
	{1./36. * M_PI, 1./9. * M_PI},
	{1./36. * M_PI, 1./9. * M_PI},
	/* long rectangle, no meridian */
	{35./18. * M_PI, -1./18. * M_PI},
	{35./18. * M_PI, -1./18. * M_PI},
	{35./18. * M_PI, -1./18. * M_PI},
	/* long rectangle, meridian */
	{17./18. * M_PI, 1./3. * M_PI},
	{17./18. * M_PI, 1./3. * M_PI},
	{17./18. * M_PI, 1./3. * M_PI},
	/* rectangle around */
	{2 * M_PI, 1/18. * M_PI},
	{2 * M_PI, 1/18. * M_PI},
	/* triangle at north pole, over meridian */
	{1./9. * M_PI, 1./2. * M_PI},
	{1./9. * M_PI, 1./2. * M_PI}};
    float pointList[NTESTS][2] = {/* small rectangle, no meridian, point in */
	{23./30. * M_PI, 13./45. * M_PI},
	/* small rectangle, no meridian, point out (right) */
	{15./18. * M_PI, 13./45. * M_PI},
	/* small rectangle, meridian, point in left */
	{-1./36. * M_PI, 0},
	/* small rectangle, meridian, point in right */
	{1./180. * M_PI, -1./90. * M_PI},
	/* small rectangle, meridian, point out (left) */
	{-1./9. * M_PI, 0},
	/* long rectangle, no meridian point in left */
	{5./6. * M_PI, -1./9. * M_PI},
	/* long rectangle, no meridian point in right */
	{5./3. * M_PI, -1./9. * M_PI},
	/* long rectangle, no meridian point out (center, top) */
	{0, 0},
	/* long rectangle, meridian, point in left */
	{-16./18. * M_PI, 25./180. * M_PI},
	/* long rectangle, meridian, point in right */
	{16./18. * M_PI, 25./180. * M_PI},
	/* long rectangle, meridian, point out (bottom) */
	{1./2. * M_PI, 0},
	/* rectangle around point in*/
	{0, 0},
	/* rectangle around, point out (top) */
	{1, 1},
	/* triangle at north pole, over meridian, point in */
	{0, 8/18. * M_PI},
	/* triangle at north pole, over meridian, point out (right) */
	{2./9. * M_PI, 8./18. * M_PI}};
    short int truthList[NTESTS] = {
	/* small rectangle, no meridian, point out */
	1, 
	/* small rectangle, no meridian, point out */
	0,
	/* small rectangle, meridian, point in left */
	1,
	/* small rectangle, meridian, point in right */
	1,
	/* small rectangle, meridian, point out */
	0,
	/* long rectangle, no meridian point in left */
	1,
	/* long rectangle, no meridian point in right */
	1,
	/* long rectangle, no meridian point out (center, top) */
	0,
	/* long rectangle, meridian, point in left */
	1,
	/* long rectangle, meridian, point in right */
	1,
	/* long rectangle, meridian, point out (bottom) */
	0,
	/* rectangle around point in*/
	1,
	/* rectangle around point out (top) */
	0,
	/* triangle at north pole, over meridian, point in */
	1,
	/* triangle at north pole, over meridian, point out (right) */
	0};

    int i;
    kd_point *point;
    kdNode *node;

    point = kd_malloc(sizeof(kd_point), "kd_point");
    point->point = kd_malloc(2 * sizeof(float), "point position");
    point->data = NULL;
    for(i = 0; i < NTESTS; i++) {
	point->point[0] = pointList[i][0];
	point->point[1] = pointList[i][1];
	node = kd_allocNode(point, 0, point->point, point->point, 0, 2);
	if (truthList[i] != kd_sph_isPointInRect(node, minList[i], 
						 maxList[i])) {
	    fprintf(stderr, "Test %d, Point (%.3f, %.3f)\n", i, 
		    point->point[0], point->point[1]);
	    fprintf(stderr, 
		    "Rectangle corners min (%.3f, %.3f), max (%.3f, %.3f)\n",
		    minList[i][0], minList[i][1], maxList[i][0], 
		    maxList[i][1]);
	    fprintf(stderr, "Expected truth value: %d\n", truthList[i]);
	    exit(EXIT_FAILURE);
	}
	kd_destroyTree(node, NULL);
    }
    return 0;
}

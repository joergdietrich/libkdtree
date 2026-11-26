#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <kdtree.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define NPOINTS 100000

int
main()
{
    float p[2], q[2];
    float d;
    int i;

    /* The distance of a point to itself is exactly zero */
    srand48(11);
    for(i = 0; i < NPOINTS; i++) {
        p[0] = drand48() * 2 * M_PI;
        p[1] = asin(2 * drand48() - 1);
        if ((d = kd_sph_dist_sq(p, p)) != 0) {
            fprintf(stderr, "Point (%.7f, %.7f) has distance %g to itself\n",
                    p[0], p[1], d);
            exit(EXIT_FAILURE);
        }
    }

    /* Pole to equator */
    p[0] = 1;
    p[1] = M_PI / 2;
    q[0] = 2;
    q[1] = 0;
    if (fabs(kd_sph_dist(p, q) - M_PI / 2) > 1e-6)
        exit(EXIT_FAILURE);

    /* Antipodal points */
    p[0] = 0.5;
    p[1] = 0.3;
    q[0] = 0.5 + M_PI;
    q[1] = -0.3;
    if (fabs(kd_sph_dist(p, q) - M_PI) > 1e-3)
        exit(EXIT_FAILURE);

    /* Small separation along a meridian */
    p[0] = q[0] = 1;
    p[1] = 0.3;
    q[1] = 0.3 + 1e-5;
    if (fabs(kd_sph_dist(p, q) / (q[1] - p[1]) - 1) > 1e-3)
        exit(EXIT_FAILURE);

    return 0;
}

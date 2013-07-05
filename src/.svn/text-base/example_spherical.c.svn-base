#include "../config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TIME_WITH_SYS_TIME
#include <sys/time.h>
#endif
#include <time.h>

#include <kdtree.h>


static double sqrarg __attribute__ ((unused));
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#ifdef HAVE_GETTIMEOFDAY
unsigned int
get_microsec(void)
{
    struct timeval timeval;

    gettimeofday(&timeval, 0);

    return timeval.tv_sec * 1000000 + timeval.tv_usec;
}
#else
unsigned int
get_microsec(void)
{

    return time(NULL) * 1000000;
}
#endif

void
printData(void *data)
{
    unsigned long *a = (unsigned long *) data;
    printf("Node data: %ld\n", *a);
}


void *
data_constr(void *data)
{
    unsigned long *a;

    a = malloc(sizeof(unsigned long));
    memcpy(a, (unsigned long *) data, sizeof(unsigned long));
    return a;
}

void
data_destr(void *data)
{
    free(data);
}

int
main(int argc, char **argv)
{
    unsigned long i;
    unsigned int npoints;
    int nthreads;
    int j;
    unsigned int start, buildTime, searchTime, rectSearchTime;
    struct kd_point *pointlist;
    struct kdNode *kdTree;
    struct kdNode *nearest;
    float rect_min[2], rect_max[2];
    float min[2], max[2];
    float range;
    FILE *f;
    struct pqueue *result;
    struct resItem *p;

    if (argc != 3) {
        fprintf(stderr,
                "Usage: example_cartesian <no of points> <no of threads>\n");
        exit(EXIT_FAILURE);
    }
    npoints = atoi(argv[1]);
    nthreads = atoi(argv[2]);
    srand48(11);

    /*
     * The kd-tree is build from an array of kd_points 
     */
    pointlist = malloc(npoints * sizeof(struct kd_point));

    /*
     * Cover the entire sphere with points
     */
    min[0] = 0;
    min[1] = -M_PI / 2.;
    max[0] = 2 * M_PI;
    max[1] = M_PI / 2.;

#ifndef HAVE_GETTIMEOFDAY
    fprintf(stderr,
            "Your system does not support the gettimeofday() function.\n"
            "The time resolution of the available time() function may be\n"
            "too low for meaningful results. The library will still work\n"
            "without problems.\n\n");
#endif

    /*
     * We randomly create points and data and write the list of points
     * to a file 
     */
    printf("Creating random points ...");
    fflush(stdout);
    f = fopen("pointlist_test.dat", "w");
    for(i = 0; i < npoints; i++) {
        pointlist[i].point = malloc(2 * sizeof(float));
        /*
         * In this case the data container of the kd_points is a
         * pointer to an unsigned long simply containing the running
         * number of this point. 
         */
        pointlist[i].data = (unsigned long *) malloc(sizeof(unsigned long));
        memcpy(pointlist[i].data, &i, sizeof(unsigned long));
        for(j = 0; j < 2; j++) {
            pointlist[i].point[j] = drand48() * (max[j] - min[j]) + min[j];
            if (j == 0 && pointlist[i].point[j] < 0)
                pointlist[i].point[j] += 2 * M_PI;
        }
        fprintf(f, "%f\t%f\t%ld\n", pointlist[i].point[0],
                pointlist[i].point[1], *(unsigned long *) pointlist[i].data);
    }
    fclose(f);
    printf("done.\nYou can find the point list in pointlist_test.dat\n\n");

    /*
     * Now we build the kd-tree from the pointlist 
     */
    printf("Building the kd-tree ... ");
    fflush(stdout);
    start = get_microsec();
    if ((kdTree =
         kd_sph_buildTree(pointlist, npoints, data_constr, data_destr, min,
                          max, nthreads)) == NULL) {
        fprintf(stderr, "Error building kd-tree\n");
        exit(EXIT_FAILURE);
    }
    buildTime = get_microsec() - start;
    printf("done (in %.3e seconds).\n\n", buildTime / 1e6);

    /*
     * Now that the tree is built we can free the kd_point array 
     */
    for(i = 0; i < npoints; i++) {
        free(pointlist[i].data);
        free(pointlist[i].point);
    }
    free(pointlist);


    /*
     * Now we perform an orthogonal range search. The results will be
     * in a min-max heap. There are special functions to retrieve
     * points from the heap if the order is important (e.g. in range
     * and distance searches). Here the order is unspecified and it is
     * faster to simply loop over the heap. 
     */

    /*
     * Define the rectangle corners 
     */
    rect_min[0] = -1.15;
    rect_min[1] = -0.5;
    rect_max[0] = 1.;
    rect_max[1] = 0.5;

    printf("Searching for points in the rectangle defined by corners\n"
           "(%.2f, %.2f), (%.2f, %.2f) ... ", rect_min[0], rect_min[1],
           rect_max[0], rect_max[1]);
    fflush(stdout);
    start = get_microsec();

    if ((result = kd_sph_ortRangeSearch(kdTree, rect_min, rect_max)) == NULL) {
        fprintf(stderr, "Orthogonal range search failed.\n");
        exit(EXIT_FAILURE);
    }

    rectSearchTime = get_microsec() - start;
    printf("done (in %.3e seconds)\n", rectSearchTime / 1e6);
    printf("Found %d points in reactangle.\n", result->size - 1);

    /*
     * Write the results to a file 
     */
    f = fopen("ort_range_result.dat", "w");
    /*
     * Loop over the heap. Note that the element with index 0 is
     * unused. We start at 1. 
     */
    for(i = 1; i < result->size; i++) {
        fprintf(f, "%f\t%f\n", result->d[i]->node->location[0],
                result->d[i]->node->location[1]);
        /*
         * Do not forget to free the memory. First free the heap
         * element. Note that at this point the heap is destroyed. It
         * does not fulfill the heap property anymore. 
         */
        free(result->d[i]);
    }
    fclose(f);
    /*
     * Free the heap 
     */
    free(result->d);
    /*
     * Free the heap information structure 
     */
    free(result);

    printf("You can find the point list in ort_range_result.dat\n\n");

    /*
     * Now we do a nearest neighbor search 
     */
    range = SQR(2 * M_PI);      /* This has to be bigger than the presumed
                                 * maximum distance to the NN but smaller
                                 * than once around the sphere. The content
                                 * of this variable is replaced with the
                                 * distance to the NN squared. */
    rect_min[0] = 1.5;
    rect_min[1] = 0.12920;
    printf
        ("Searching for nearest neighbor around the point (%.2f, %.2f) ... ",
         rect_min[0], rect_min[1]);
    fflush(stdout);
    start = get_microsec();

    nearest = kd_sph_nearest(kdTree, rect_min, &range);

    searchTime = get_microsec() - start;
    printf("\n done (in %.3e seconds).\nDistance to the nearest neigbor is "
           "%.5e.\n"
           "This is the node containing the neareast neighbor:\n",
           searchTime / 1e6, sqrt(range));
    kd_printNode(nearest, printData);


    /*
     * Perform a range search around a point. Return the results
     * ordered by distance. The results are in a priority
     * queue. Retrieving each result from the queue takes O(log n)
     * time, where n is the size of the queue. 
     * 
     * If the results are unordered the cost of building the queue is
     * O(1) and one could simply loop over it (with O(n) cost), like
     * it is done for the orthogonal range search above, instead of
     * using the queue removal (which costs O(n log n)).
     */
    range = SQR(1.50);          /* range contains the search radius squared */

    /*
     * Search around the zero meridian 
     */
    rect_max[0] = 0;
    rect_max[1] = 0.3;

    printf("Performing range search with radius %.2f around point "
           "(%.2f, %.2f) ... ", sqrt(range), rect_max[0], rect_max[1]);
    fflush(stdout);
    start = get_microsec();

    if ((result = kd_sph_range(kdTree, rect_max, &range, KD_ORDERED)) == NULL) {
        fprintf(stderr, "Range search failed.\n");
        exit(EXIT_FAILURE);
    }

    searchTime = get_microsec() - start;
    j = result->size - 1;
    printf("\n done (in %.3e seconds). Found %d points.\n"
           "Results are ordered by distance (ascending order).\n",
           searchTime / 1e6, j);
    f = fopen("range_result.dat", "w");
    /*
     * pqremove_min returns the points in ascending order. Use
     * pqremove_max to get the points in decending order. 
     */
    while (pqremove_min(result, &p)) {
        fprintf(f, "%f\t%f\t%g\n", p->node->location[0], p->node->location[1],
                p->dist_sq);
        /*
         * Free the result node taken from the heap 
         */
        free(p);
    }
    fclose(f);
    /*
     * free the unused first heap element 
     */
    free(result->d);
    /*
     * and free the heap information structure 
     */
    free(result);
    printf("You can find the point list in range_result.dat.\n\n");

    /*
     * Now we do a q-nearest neighbor search 
     */
    range = SQR(2);             /* do not forget to reset the range. It was
                                 * modified by kd_range above. */
    printf("Searching for the %d nearest neighbors around the point "
           "(%.2f, %.2f) ... ", j, rect_max[0], rect_max[1]);
    fflush(stdout);
    start = get_microsec();

    if ((result = kd_sph_qnearest(kdTree, rect_max, &range, j)) == NULL) {
        fprintf(stderr, "Q-nearest neighbor search failed.\n");
        exit(EXIT_FAILURE);
    }

    searchTime = get_microsec() - start;
    printf("\n done (in %.3e seconds). Found %d points.\n"
           "Results are ordered by distance (descending order).\n",
           searchTime / 1e6, result->size - 1);

    f = fopen("qnearest_result.dat", "w");
    /*
     * pqremove_min returns the points in ascending order. Use
     * pqremove_max to get the points in decending order. 
     */
    while (pqremove_max(result, &p)) {
        fprintf(f, "%f\t%f\t%g\n", p->node->location[0], p->node->location[1],
                p->dist_sq);
        /*
         * Free the result node taken from the heap 
         */
        free(p);
    }
    fclose(f);
    /*
     * free the heap 
     */
    free(result->d);
    /*
     * and free the heap information structure 
     */
    free(result);
    printf("You can find the point list in qnearest_result.dat\n"
           "(these should be the exact same points as in range_result.dat\n"
           "but with order reversed).\n\n");

    /*
     * Free memory occupied by kd-tree 
     */
    kd_destroyTree(kdTree, data_destr);
    return 0;
}

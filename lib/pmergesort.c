#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct param_t {
    char *a;
    char *b;
    size_t first;
    size_t nmemb;
    size_t size;
    int (*cmp) (const void *, const void *, int axis);
    int axis;
    int max_threads;
} param_t;

extern void qsortR(const void *base0, size_t n, size_t size,
                   int (*compar) (const void *, const void *, int axis),
                   int axis);


void pm_buildparams(struct param_t *p, void *a, void *b, size_t first,
                    size_t nmemb, size_t size,
                    int (*cmp) (const void *, const void *, int),
                    int axis, int max_threads);
int pmergesort(void *base, size_t nmemb, size_t size,
               int (*compar) (const void *, const void *, int),
               int axis, int max_threads);
void *mergesort_t(void *args);


int
pmergesort(void *base, size_t nmemb, size_t size,
           int (*compar) (const void *, const void *, int),
           int axis, int max_threads)
{
    void *tmp;
    param_t args;

    if (nmemb < 2)
        return 1;
    if ((tmp = calloc(nmemb, size)) == NULL) {
        perror("malloc");
        return 0;
    }
    args.a = (char *) base;
    args.b = (char *) tmp;
    args.first = 0;
    args.nmemb = nmemb;
    args.size = size;
    args.cmp = compar;
    args.axis = axis;
    args.max_threads = max_threads;

    mergesort_t(&args);

    free(tmp);
    return 1;
}

void
pm_buildparams(struct param_t *p, void *a, void *b, size_t first,
               size_t nmemb, size_t size,
               int (*cmp) (const void *, const void *, int),
               int axis, int max_threads)
{

    p->a = a;
    p->b = b;
    p->first = first;
    p->nmemb = nmemb;
    p->size = size;
    p->cmp = cmp;
    p->axis = axis;
    p->max_threads = max_threads;
}


void *
mergesort_t(void *args)
{
    size_t i, li, ri;
    struct param_t larg, rarg;
    pthread_t thr[2];
    int threaded[2];
    param_t *mya = (param_t *) args;

    if (mya->max_threads < 2) {
        /*
         * Reached maximum number of threads allocated to this
         * branch. Proceed with sequential sort of this chunk. 
         */
        qsortR(mya->a + mya->first * mya->size, mya->nmemb, mya->size,
               mya->cmp, mya->axis);
    } else {
        /*
         * Sort each half of array a in a new thread. A half whose
         * thread cannot be started is sorted in this thread.
         */
        pm_buildparams(&larg, mya->a, mya->b, mya->first, mya->nmemb / 2,
                       mya->size, mya->cmp, mya->axis, mya->max_threads / 2);
        pm_buildparams(&rarg, mya->a, mya->b, mya->first + mya->nmemb / 2,
                       mya->nmemb - mya->nmemb / 2, mya->size,
                       mya->cmp, mya->axis, mya->max_threads / 2);

        threaded[0] = !pthread_create(&thr[0], NULL, mergesort_t, &larg);
        threaded[1] = !pthread_create(&thr[1], NULL, mergesort_t, &rarg);
        if (threaded[0])
            pthread_join(thr[0], NULL);
        else
            mergesort_t(&larg);
        if (threaded[1])
            pthread_join(thr[1], NULL);
        else
            mergesort_t(&rarg);

        /*
         * Merge the two sorted chunks of array a into array b 
         */
        li = larg.first;
        ri = rarg.first;
        for(i = mya->first; i < mya->first + mya->nmemb; i++) {
            if (li >= larg.first + larg.nmemb) {
                /*
                 * We already copied everything from the left chunk,
                 * now copy from the right 
                 */
                memcpy(mya->b + i * mya->size, mya->a + ri * mya->size,
                       mya->size);
                ri++;
            } else if (ri >= rarg.first + rarg.nmemb) {
                /*
                 * We already copied everything from the right chunk,
                 * now copy from the left 
                 */
                memcpy(mya->b + i * mya->size, mya->a + li * mya->size,
                       mya->size);
                li++;
            }
            /*
             * We can still copy from both chunks, copy the smaller
             * element 
             */
            else if (mya->cmp(mya->a + li * mya->size,
                              mya->a + ri * mya->size, mya->axis) < 1) {
                memcpy(mya->b + i * mya->size, mya->a + li * mya->size,
                       mya->size);
                li++;
            } else {
                memcpy(mya->b + i * mya->size, mya->a + ri * mya->size,
                       mya->size);
                ri++;
            }
        }
        /*
         * Now b is sorted, copy it back to a 
         */
        memcpy(mya->a + mya->size * mya->first,
               mya->b + mya->size * mya->first, mya->size * mya->nmemb);
    }
    return NULL;
}

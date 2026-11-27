/*!\file pqueue.c
 * \brief Two-ended priority queues (min-max heap)
 *
 * Implementation of a min-max heap (two-ended priority queue) as
 * introduced by Atkinson et al. (1986), Communications of the ACM 10,
 * 996.
 */

#include "kdtree_internal.h"

#define PQPRIO(p) (p->dist_sq)

static int floorLog2(uint32_t n);
static int is_max_level(int i);
static uint32_t get_max_index(struct pqueue *q);
static uint32_t pq_get_min_child_index(struct pqueue *q, uint32_t i);
static uint32_t pq_get_max_child_index(struct pqueue *q, uint32_t i);
static void pq_swap_nodes(struct pqueue *q, uint32_t i, uint32_t j);
static void trickle_down(struct pqueue *q, uint32_t i);
static void trickle_down_min(struct pqueue *q, uint32_t i);
static void trickle_down_max(struct pqueue *q, uint32_t i);
static void bubble_up(struct pqueue *q, uint32_t i);
static void bubble_up_min(struct pqueue *q, uint32_t i);
static void bubble_up_max(struct pqueue *q, uint32_t i);



/**********************************************************************
 *
 *  We start with some general helper functions 
 *
 */

/*
   Returns the floor form of a binary logarithm for a 32 bit integer.
   -1 is returned if n is 0.
   Ripped off wikipedia.
*/
static inline int
floorLog2(uint32_t n)
{
    uint32_t pos = 0;
    if (n >= 1 << 16) {
        n >>= 16;
        pos += 16;
    }
    if (n >= 1 << 8) {
        n >>= 8;
        pos += 8;
    }
    if (n >= 1 << 4) {
        n >>= 4;
        pos += 4;
    }
    if (n >= 1 << 2) {
        n >>= 2;
        pos += 2;
    }
    if (n >= 1 << 1)
        pos += 1;
    if (n == 0)
        return -1;
    return pos;
}


/* 
   Returns 1 if i is a max-sorted level, return 0 otherwise
   (min-sorted level). 
*/
static inline int
is_max_level(int i)
{
    return floorLog2(i) % 2;
}



/* Swap nodes i and j in the priority queue q */
static void
pq_swap_nodes(struct pqueue *q, uint32_t i, uint32_t j)
{
    struct resItem *tmp;

    tmp = q->d[i];
    q->d[i] = q->d[j];
    q->d[j] = tmp;
}


/* Return the array index of the maximum node */
static uint32_t
get_max_index(struct pqueue *q)
{

    uint32_t i;
    if (!q)
        return 0;

    switch (q->size) {
    case 1:{
            /*
             * Empty queue 
             */
            i = 0;
            break;
        }
    case 2:{
            /*
             * this is the only element (min and max) 
             */
            i = 1;
            break;
        }
    case 3:{
            /*
             * This is the only element on the max level 
             */
            i = 2;
            break;
        }
    default:{
            /*
             * Return the greater of the two elements on the max level 
             */
            if (PQPRIO(q->d[2]) > PQPRIO(q->d[3]))
                i = 2;
            else
                i = 3;
        }
    }
    return i;
}


/* End of general helper functions */

/* 
 *
 *  Functions to support insertion 
 *
 */

/* Move a node up the tree */
static void
bubble_up_min(struct pqueue *q, uint32_t i)
{
    /*
     * if node has a grandparent 
     */
    if (i > 3 && PQPRIO(q->d[i]) < PQPRIO(q->d[i / 4])) {
        pq_swap_nodes(q, i, i / 4);
        bubble_up_min(q, i / 4);
    }
}


static void
bubble_up_max(struct pqueue *q, uint32_t i)
{
    /*
     * if node has a grandparent 
     */
    if (i > 3 && PQPRIO(q->d[i]) > PQPRIO(q->d[i / 4])) {
        pq_swap_nodes(q, i, i / 4);
        bubble_up_max(q, i / 4);
    }
}


static void
bubble_up(struct pqueue *q, uint32_t i)
{
    if (!is_max_level(i)) {
        if (i > 1 && PQPRIO(q->d[i]) > PQPRIO(q->d[i / 2])) {
            pq_swap_nodes(q, i, i / 2);
            bubble_up_max(q, i / 2);
        } else {
            bubble_up_min(q, i);
        }
    } else {
        if (i > 1 && PQPRIO(q->d[i]) < PQPRIO(q->d[i / 2])) {
            pq_swap_nodes(q, i, i / 2);
            bubble_up_min(q, i / 2);
        } else {
            bubble_up_max(q, i);
        }
    }
}


/*
 *
 * Functions to support node removal
 *
 */

/* Get index of the smallest child or grandchild of q->d[i].
   Caller must ensure that q->d[i] has at least one child. */
static uint32_t
pq_get_min_child_index(struct pqueue *q, uint32_t i)
{
    uint32_t m;

    /*
     * First Child 
     */
    m = 2 * i;
    /*
     * Second Child 
     */
    if (q->size > 2 * i + 1)
        if (PQPRIO(q->d[m]) > PQPRIO(q->d[2 * i + 1]))
            m = 2 * i + 1;
    /*
     * First child of first child 
     */
    if (q->size > 4 * i)
        if (PQPRIO(q->d[m]) > PQPRIO(q->d[4 * i]))
            m = 4 * i;
    /*
     * Second child of first child 
     */
    if (q->size > 4 * i + 1)
        if (PQPRIO(q->d[m]) > PQPRIO(q->d[4 * i + 1]))
            m = 4 * i + 1;
    /*
     * First Child of second child 
     */
    if (q->size > 2 * (2 * i + 1))
        if (PQPRIO(q->d[m]) > PQPRIO(q->d[2 * (2 * i + 1)]))
            m = 2 * (2 * i + 1);
    /*
     * Second child of second child 
     */
    if (q->size > 2 * (2 * i + 1) + 1)
        if (PQPRIO(q->d[m]) > PQPRIO(q->d[2 * (2 * i + 1) + 1]))
            m = 2 * (2 * i + 1) + 1;
    return m;
}

/* Get index of the largest children and grandchildren of q->d[i].
   Caller must ensure that q->d[i] has at least one child. */
static uint32_t
pq_get_max_child_index(struct pqueue * q, uint32_t i)
{
    uint32_t m;

    /*
     * First Child 
     */
    m = 2 * i;
    /*
     * Second Child 
     */
    if (q->size > 2 * i + 1)
        if (PQPRIO(q->d[m]) < PQPRIO(q->d[2 * i + 1]))
            m = 2 * i + 1;
    /*
     * First child of first child 
     */
    if (q->size > 4 * i)
        if (PQPRIO(q->d[m]) < PQPRIO(q->d[4 * i]))
            m = 4 * i;
    /*
     * Second child of first child 
     */
    if (q->size > 4 * i + 1)
        if (PQPRIO(q->d[m]) < PQPRIO(q->d[4 * i + 1]))
            m = 4 * i + 1;
    /*
     * First Child of second child 
     */
    if (q->size > 2 * (2 * i + 1))
        if (PQPRIO(q->d[m]) < PQPRIO(q->d[2 * (2 * i + 1)]))
            m = 2 * (2 * i + 1);
    /*
     * Second child of second child 
     */
    if (q->size > 2 * (2 * i + 1) + 1)
        if (PQPRIO(q->d[m]) < PQPRIO(q->d[2 * (2 * i + 1) + 1]))
            m = 2 * (2 * i + 1) + 1;
    return m;
}


/* Move a node down the tree */
static void
trickle_down(struct pqueue *q, uint32_t i)
{
    if (is_max_level(i))
        trickle_down_max(q, i);
    else
        trickle_down_min(q, i);
}


static void
trickle_down_max(struct pqueue *q, uint32_t i)
{
    uint32_t m;

    /*
     * if A[i] has children 
     */
    if (q->size > 2 * i) {
        /*
         * m := Index of the largest child or grandchild of A[i] 
         */
        m = pq_get_max_child_index(q, i);
        /*
         * if A[m] is a grandchild of A[i] 
         */
        if (m / 4 == i) {
            if (PQPRIO(q->d[m]) > PQPRIO(q->d[i]))
                pq_swap_nodes(q, i, m);
            if (PQPRIO(q->d[m]) < PQPRIO(q->d[m / 2]))
                pq_swap_nodes(q, m, m / 2);
            trickle_down_max(q, m);
        } else {
            /*
             * A[m] is a child of A[i] 
             */
            if (PQPRIO(q->d[m]) > PQPRIO(q->d[i]))
                pq_swap_nodes(q, i, m);
        }
    }
}

static void
trickle_down_min(struct pqueue *q, uint32_t i)
{
    uint32_t m;

    /*
     * if A[i] has children 
     */
    if (q->size > 2 * i) {
        /*
         * m := Index of the smallest child or grandchild of A[i] 
         */
        m = pq_get_min_child_index(q, i);
        /*
         * if A[m] is a grandchild of A[i] 
         */
        if (m / 4 == i) {
            if (PQPRIO(q->d[m]) < PQPRIO(q->d[i]))
                pq_swap_nodes(q, i, m);
            if (PQPRIO(q->d[m]) > PQPRIO(q->d[m / 2]))
                pq_swap_nodes(q, m, m / 2);
            trickle_down_min(q, m);
        } else {
            /*
             * A[m] is a child of A[i] 
             */
            if (PQPRIO(q->d[m]) < PQPRIO(q->d[i]))
                pq_swap_nodes(q, i, m);
        }
    }
}




/* 
 ***********************************************************************
 *
 *                USER FUNCTIONS
 *
 ***********************************************************************
 */

/*
 * Initializes the priority queue q, or allocates a new one if q is
 * NULL. Memory for n items is allocated up front, and the queue grows
 * by another n items whenever it fills up. An n of 0 counts as 1.
 * Returns the queue, or NULL on error.
 */
struct pqueue *
pqinit(struct pqueue *q, uint32_t n)
{
    struct pqueue *tmp = q;

    /* The queue grows in steps of n, so n must not be 0 */
    if (n == 0)
        n = 1;
    if (!q && !(q = malloc(sizeof(struct pqueue)))) {
        return NULL;
    }
    if (!(q->d = malloc(sizeof(struct resItem *) * n))) {
        if (!tmp)
            free(q);
        return NULL;
    }
    q->avail = q->step = n;
    q->size = 1;
    return q;
}


/*
 * Inserts the item d into the queue q. Returns 1 on success, or 0 if
 * q is NULL or the queue could not grow.
 */
int
pqinsert(struct pqueue *q, struct resItem *d)
{
    struct resItem **tmp;
    uint32_t i, newsize;

    if (!q)
        return 0;

    /*
     * allocate more memory if necessary 
     */
    if (q->size >= q->avail) {
        newsize = q->size + q->step;
        if (!(tmp = realloc(q->d, sizeof(struct resItem *) * newsize))) {
            return 0;
        };
        q->d = tmp;
        q->avail = newsize;
    }

    /*
     * insert item 
     */
    i = q->size++;
    q->d[i] = d;
    bubble_up(q, i);
    return 1;
}


/*!
 * \brief remove the highest-ranking (minimum) item from the queue.
 *
 * \param q a pointer to a priority queue.
 *  
 * \param d a pointer to the struct resItem * variable that will hold the datum
 * corresponding to the queue item removed.
 *
 * \return non-NULL if an item has been removed. The variable that d
 * points to now contains the datum associated with the item in
 * question; or NULL if item could be removed. Either the queue
 * pointer provided was NULL, or the queue was empty. The chunk of
 * memory that d points to has not been modified.
 */
struct resItem **
pqremove_min(struct pqueue *q, struct resItem **d)
{
    if (!q || q->size == 1)
        return NULL;
    *d = q->d[1];
    q->d[1] = q->d[--q->size];
    trickle_down(q, 1);
    return d;
}


/*!
 * \brief remove the lowest-ranking (maximum) item from the queue.
 *
 * \param q a pointer to a priority queue.
 *
 * \param d a pointer to the struct resItem * variable that will hold the datum
 * corresponding to the queue item removed.
 *
 * \return non-NULL if an item has been removed. The variable that d
 * points to now contains the datum associated with the item in
 * question; or NULL if item could be removed. Either the queue
 * pointer provided was NULL, or the queue was empty. The chunk of
 * memory that d points to has not been modified.
 */
struct resItem **
pqremove_max(struct pqueue *q, struct resItem **d)
{
    uint32_t i;

    if (!q || q->size == 1)
        return NULL;
    i = get_max_index(q);
    *d = q->d[i];
    q->d[i] = q->d[--q->size];
    trickle_down(q, i);
    return d;
}


/*!
 * \brief access highest-ranking (minimum) item without removing it.
 *
 * \param q a pointer to a priority queue.
 *
 * \param d a pointer to the struct resItem * variable that will hold the datum
 * corresponding to the highest-ranking item.
 *                
 * \return non-NULL in case of success. The variable that d points to
 * now contains the datum associated with the highest-ranking item;
 * NULL in case of failure. Either the queue pointer provided was
 * NULL, or the queue was empty. The chunk of memory that d points to
 * has not been modified.
 */
struct resItem **
pqpeek_min(struct pqueue *q, struct resItem **d)
{
    if (!q || q->size == 1)
        return NULL;
    *d = q->d[1];
    return d;
}


/*!
 * \brief access lowest-ranking (maximum) item without removing it.
 *
 * \param q a pointer to a priority queue.
 *
 * \param d a pointer to the struct resItem * variable that will hold the datum
 * corresponding to the highest-ranking item.
 *                
 * \return non-NULL in case of success. The variable that d points to
 * now contains the datum associated with the highest-ranking item;
 * NULL in case of failure. Either the queue pointer provided was
 * NULL, or the queue was empty. The chunk of memory that d points to
 * has not been modified.
 */
struct resItem **
pqpeek_max(struct pqueue *q, struct resItem **d)
{
    if (!q || q->size == 1)
        return NULL;
    *d = q->d[get_max_index(q)];
    return d;
}


/*!
 * \brief free a result queue and the items still in it
 *
 * \param q a queue returned by a search, or NULL. Items already taken
 * out with pqremove_min() or pqremove_max() are the caller's to free.
 */
void
pqfree(struct pqueue *q)
{
    uint32_t i;

    if (!q)
        return;
    for(i = 1; i < q->size; i++)
        free(q->d[i]);
    free(q->d);
    free(q);
}

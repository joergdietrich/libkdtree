#include "kdtree.h"

#ifndef _PQUEUE_H
#define _PQUEUE_H

#include <stdint.h>
#include <stdlib.h>

#define PQPRIO(p) (p->dist_sq)

int floorLog2(uint32_t n);
int is_max_level(int i);
uint32_t get_max_index(struct pqueue *q);
uint32_t pq_get_min_child_index(struct pqueue *q, uint32_t i);
uint32_t pq_get_max_child_index(struct pqueue *q, uint32_t i);
void pq_swap_nodes(struct pqueue *q, uint32_t i, uint32_t j);
void trickle_down(struct pqueue *q, uint32_t i);
void trickle_down_min(struct pqueue *q, uint32_t i);
void trickle_down_max(struct pqueue *q, uint32_t i);
void bubble_up(struct pqueue *q, uint32_t i);
void bubble_up_min(struct pqueue *q, uint32_t i);
void bubble_up_max(struct pqueue *q, uint32_t i);

#endif

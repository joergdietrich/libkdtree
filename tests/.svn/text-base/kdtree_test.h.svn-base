struct kd_point *create_pointlist(unsigned long npoints, float *min, 
				  float *max, unsigned short dim);
int points_eq(struct kdNode *nkd, float *nn, unsigned short dim);
int sorted_queues_eq(struct pqueue *q1, struct pqueue *q2, unsigned short dim);
int  unsorted_queues_eq(struct pqueue *p1, struct pqueue *p2, 
			unsigned short dim);
int  is_pointInRect(struct kd_point point, float *rect_min, float *rect_max,
	       unsigned short dim);
struct pqueue *naive_ortRangeSearch(struct kd_point *pointlist, 
				    unsigned int npoints, float *rect_min, 
				    float *rect_max, unsigned short dim);
int sph_is_pointInRect(struct kd_point point, float *rect_min, float *rect_max,
		       unsigned short dim);
struct pqueue *naive_sph_ortRangeSearch(struct kd_point *pointlist, 
					unsigned int npoints, float *rect_min,
					float *rect_max, unsigned short dim);

#define EPS 1e-7


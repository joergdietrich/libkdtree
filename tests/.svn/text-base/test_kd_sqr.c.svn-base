#include <stdlib.h>
#include <kdtree.h>

int 
main()
{
    float in[3] = {0, 2, -2};
    float out[3] = {0, 4, 4};
    int i;

    for(i = 0; i < 3; i++) {
	if (kd_sqr(in[i]) != out[i])
	    exit(EXIT_FAILURE);
    }
    return 0;
}

    

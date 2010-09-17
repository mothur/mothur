#ifndef GUARD_fisher2
#define GUARD_fisher2

#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <limits.h> 


#ifdef __cplusplus
extern "C" {
#endif

#define SINT_MAX INT_MAX

#define	max(a, b)		((a) < (b) ? (b) : (a))
#define	min(a, b)		((a) > (b) ? (b) : (a))


void fexact(int *, int *, double *, int *,
       double *, double *, double *, double *,
       double *, int *);
	   
#ifdef __cplusplus	   
}
#endif

#endif 


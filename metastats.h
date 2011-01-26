#ifndef METASTATS2
#define METASTATS2

/*
 *  metastats.h
 *  Mothur
 *
 *  Created by westcott on 9/16/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "fisher2.h"

void testp(double *permuted_ttests,int *B,double *permuted,double 
	   *Imatrix,int *nc,int *nr,int *g,double *Tinitial,double *ps);
void permute_matrix(double *Imatrix,int *nc,int *nr,double 
		    *permuted,int *g,double *trial_ts,double *Tinitial,double 
		    *counter);
void permute_array(int *array, int n);
void calc_twosample_ts(double *Pmatrix,int *g,int *nc,int *nr,double 
		       *Ts,double *Tinitial,double *counter1);
void meanvar(double *pmatrix,int *g,int *nr,int *nc,double *storage);
void start(double *Imatrix,int *g,int *nr,int *nc,double *testing,
			 double storage[][9]);

int metastat_main (char*, int, int, double, int, double**, int);

#ifdef __cplusplus	   
}
#endif

#endif




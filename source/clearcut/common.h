/*
 * common.h
 *
 * $Id$
 *
 *****************************************************************************
 *
 * Copyright (c) 2004,  Luke Sheneman
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions 
 * are met:
 * 
 *  + Redistributions of source code must retain the above copyright 
 *    notice, this list of conditions and the following disclaimer. 
 *  + Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in 
 *    the documentation and/or other materials provided with the 
 *    distribution. 
 *  + The names of its contributors may not be used to endorse or promote 
 *    products derived  from this software without specific prior 
 *    written permission. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.  
 *
 *****************************************************************************
 *
 * A header file filled with common definitions and simple inline functions
 *
 *****************************************************************************
 *
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_NJ_COMMON_H_
#define _INC_NJ_COMMON_H_ 1

#include <math.h>
#include <float.h>


#define NJ_AMBIGUITY_CHAR    63  /* ? character */


/*
 * this macro defines the number of cells in the diagonal matrix 
 * based on the number of taxa involved
 *
 */
#define NJ_NCELLS(a)       ( ((a)*(a+1))/2 )




/*
 * NJ_MAP() - 
 *
 * Thus function maps i, j coordinates to the correct offset into 
 * the distance matrix
 *
 */
static inline
long int 
NJ_MAP(long int i,
       long int j,
       long int ntaxa) {
  
  return((i*(2*ntaxa-i-1))/2 + j);
}


static inline
int
NJ_FLT_EQ(float x,
	  float y) {
  
  if(fabs(x - y)<FLT_EPSILON) {
    return(1);
  } else {
    return(0);
  }
}



static inline
int
NJ_FLT_LT(float x,
	  float y) {
  
  if(NJ_FLT_EQ(x, y)) {
    return(0);
  } else {
    if(x < y) {
      return(1);
    } else {
      return(0);
    }
  }
}


static inline
int
NJ_FLT_GT(float x,
	  float y) {
  
  if(NJ_FLT_EQ(x, y)) {
    return(0);
  } else {
    if(x > y) {
      return(1);
    } else {
      return(0);
    }
  }
}




#endif /* _INC_NJ_COMMON_H_ */




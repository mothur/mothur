/*
 * dmat.h
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
 * Distance matrix parser header file
 *
 *****************************************************************************
 *
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 */


#ifndef _INC_DMAT_H_
#define _INC_DMAT_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

#include "clearcut.h"


#define NJ_INITIAL_BUFSIZE 32

#define NJ_NAME_STATE  100
#define NJ_FLOAT_STATE 101
#define NJ_WS_STATE    102
#define NJ_EOF_STATE   103

#define NJ_PARSE_SYMMETRIC 100
#define NJ_PARSE_LOWER     101
#define NJ_PARSE_UPPER     102
#define NJ_PARSE_UNKNOWN   103


/* some data structures */
typedef struct _NJ_DIST_TOKEN_STRUCT {
  
  char *buf;
  long int bufsize;
  int type;

} NJ_DIST_TOKEN;



/* some function prototypes */

DMAT *
NJ_parse_distance_matrix(NJ_ARGS *nj_args);

void
NJ_output_matrix(NJ_ARGS *nj_args,
		 DMAT *dmat);

#ifdef __cplusplus
}
#endif

#endif /* _INC_DMAT_H_ */


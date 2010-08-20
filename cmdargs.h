/*
 * njdist.h
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
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_NJ_CMDARGS_H_
#define _INC_NJ_CMDARGS_H_ 1

#include "clearcut.h"


/* some datatypes */
typedef struct _STRUCT_NJ_ARGS {

  char *infilename;   /* the name of the input file                  */
  char *outfilename;  /* the name of the output tree                 */
  char *matrixout;    /* the name of the distance matrix output file */

  int input_mode;  
  int aligned_flag;
  
  int verbose_flag;
  int quiet_flag;

  int stdin_flag;
  int stdout_flag;
  
  int help;
  int version;

  int norandom;
  int shuffle;
  
  int dna_flag;
  int protein_flag;
  
  int seed;

  /* correction models for distance */
  int correction_model;
  int jukes_flag;
  int kimura_flag;
  
  /* flag for using traditional neighbor-joining */
  int neighbor;
  
  /* number of trees to output */
  int ntrees;
  
  /* exponential notation output */
  int expblen;  /* exp notation for tree branch lengths */
  int expdist;  /* exp notation for distances in matrix output */
  
} NJ_ARGS;



/* some function prototypes */

NJ_ARGS *
NJ_handle_args(int argc,
	       char *argv[]);

void
NJ_print_args(NJ_ARGS *nj_args);

void
NJ_usage(void);


#endif /* _INC_NJ_CMDARGS_H_ */


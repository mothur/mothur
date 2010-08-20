/*
 * fasta.h
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


#ifndef _INC_NJ_FASTA_H_
#define _INC_NJ_FASTA_H_ 1

#include "clearcut.h"

#define NJ_INITIAL_BUFFER_SIZE   512
#define NJ_INITIAL_NSEQS         64

#define NJ_FASTA_MODE_TITLE      100
#define NJ_FASTA_MODE_SEQUENCE   101
#define NJ_FASTA_MODE_NEWLINE    102
#define NJ_FASTA_MODE_UNKNOWN    103


typedef struct _STRUCT_NJ_ALIGNMENT {

  long int nseq;
  long int length;
  
  char **titles;
  
  char *data;

} NJ_alignment;


NJ_alignment *
NJ_read_fasta(NJ_ARGS *nj_args);

void
NJ_print_alignment(NJ_alignment *alignment);

void
NJ_free_alignment(NJ_alignment *alignment);

int
NJ_taxaname_unique(NJ_alignment *alignment);

#endif /* _INC_NJ_FASTA_H_ */





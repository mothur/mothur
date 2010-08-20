/*
 * dist.h
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
 * Compute a distance matrix given a set of sequences
 *
 *****************************************************************************
 *
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */


#ifndef _INC_DIST_H_
#define _INC_DIST_H_ 1

#include "fasta.h"
#include "clearcut.h"



/* 
 * An arbitrarily large distance to represent distances
 * which are too great to accurately correct.
 */
#define NJ_BIGDIST 10.0  



/* some function prototypes */
DMAT *
NJ_build_distance_matrix(NJ_ARGS *nj_args);

DMAT *
NJ_compute_dmat(NJ_ARGS *nj_args,
		NJ_alignment *alignment);


float
NJ_pw_percentid(NJ_alignment *alignment,
		long int x,
		long int y);

long int
NJ_pw_differences(NJ_alignment *alignment,
		  long int x,
		  long int y,
		  long int *residues);

void
NJ_no_correction(DMAT *dmat,
		 NJ_alignment *alignment);

void
NJ_DNA_jc_correction(DMAT *dmat,
		     NJ_alignment *alignment);

void
NJ_PROTEIN_jc_correction(DMAT *dmat,
			 NJ_alignment *alignment);

void
NJ_DNA_k2p_correction(DMAT *dmat,
		      NJ_alignment *alignment);

void
NJ_PROTEIN_kimura_correction(DMAT *dmat,
			     NJ_alignment *alignment);

void
NJ_DNA_count_tt(NJ_alignment *alignment,
		long int x,
		long int y,
		long int *transitions,
		long int *transversions,
		long int *residues);


#endif /* _INC_DIST_H_ */










/*
 * dist.c
 *
 * $Id$
 *
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "dayhoff.h"
#include "fasta.h"
#include "distclearcut.h"




/*
 * NJ_build_distance_matrix() - 
 *
 * Given a filename for an alignment, read the alignment
 * into memory and then compute the distance matrix
 * using the appropriate correction model
 */
DMAT *
NJ_build_distance_matrix(NJ_ARGS *nj_args) {
  
  DMAT *dmat;
  NJ_alignment *alignment;

  /* Read an alignment in FASTA format */
  alignment = 
    NJ_read_fasta(nj_args);
 
  if(!alignment) {
    return(nullptr);
  }

  /* 
   * Given a global multiple sequence alignment (MSA) and
   * a specified distance correction model, compute a 
   * corrected distance matrix
   *
   * From proteins, we may want to allow users to specify
   * a substitution matrix (feature)
   */

  dmat = 
    NJ_compute_dmat(nj_args,
		    alignment);

  // NJ_print_taxanames(dmat);

  if(!dmat) {
    fprintf(stderr, "Clearcut: Error computing distance matrix\n");
  }
 
  /* now free the memory associated with the alignment */
  NJ_free_alignment(alignment);

  return(dmat);
}





/* 
 * NJ_compute_dmat() - 
 *
 * Given an alignment and a correction model, compute the 
 * distance matrix and return it
 *
 */
DMAT *
NJ_compute_dmat(NJ_ARGS *nj_args,
		NJ_alignment *alignment) {

  DMAT *dmat;
  long int i;
  
  
  /* allocate distance matrix here */
  dmat = (DMAT *)calloc(1, sizeof(DMAT));
  if(!dmat) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_compute_dmat()\n");
    return(nullptr);
  }
  
  dmat->ntaxa = alignment->nseq;
  dmat->size  = alignment->nseq;

  /* allocate memory to hold the taxa names */
  dmat->taxaname = (char **)calloc(alignment->nseq, sizeof(char *));
  if(!dmat->taxaname) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_compute_dmat()\n");
    return(nullptr);
  }
  
  /* copy sequence titles */
  for(i=0;i<alignment->nseq;i++) {
    dmat->taxaname[i] = (char *)calloc(strlen(alignment->titles[i])+1, sizeof(char));
    if(!dmat->taxaname[i]) {
      fprintf(stderr, "Clearcut: Memory allocation error in NJ_compute_dmat()\n");
      return(nullptr);
    }
      
      *dmat->taxaname[i] = '\0'; strncat(dmat->taxaname[i], alignment->titles[i], strlen(alignment->titles[i])+1);
      
    //strncpy(dmat->taxaname[i], alignment->titles[i], sizeof dmat->taxaname[i] - strlen (dmat->taxaname[i]) - 1);
      
  }

  /* allocate val matrix in dmat */
  dmat->val = (float *)calloc(dmat->ntaxa*dmat->ntaxa, sizeof(float));

  if(!dmat->val) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_compute_dmat()\n");
    return(nullptr);
  }

  /* now lets allocate space for the r and r2 columns */
  dmat->r  = (float *)calloc(dmat->ntaxa, sizeof(float));
  dmat->r2 = (float *)calloc(dmat->ntaxa, sizeof(float));

  /* track some memory addresses */
  dmat->rhandle   = dmat->r;
  dmat->r2handle  = dmat->r2;
  dmat->valhandle = dmat->val;
  
  /* apply model correction to matrix */
  switch(nj_args->correction_model) {

  case NJ_MODEL_JUKES:
    
    if(nj_args->dna_flag) {
      NJ_DNA_jc_correction(dmat, alignment);
    } else if(nj_args->protein_flag) {
      NJ_PROTEIN_jc_correction(dmat, alignment);
    } else {
      fprintf(stderr, "Clearcut: Need to know sequence type for Jukes-Cantor model correction.\n");
      return(nullptr);
    }

    break;

  case NJ_MODEL_KIMURA:

    if(nj_args->dna_flag) {
      NJ_DNA_k2p_correction(dmat, alignment);
    } else if(nj_args->protein_flag) {
      NJ_PROTEIN_kimura_correction(dmat, alignment);
    } else {
      fprintf(stderr, "Clearcut: Need to know sequence type for Kimura model correction.\n");
      return(nullptr);
    }

    break;

  case NJ_MODEL_NONE:

    NJ_no_correction(dmat, alignment);

    break;

  default:
    fprintf(stderr, "Clearcut: Invalid distance correction model.\n");
    return(nullptr);
  }
 
  return(dmat);
}





/*
 * NJ_no_correction() - 
 *
 * Compute the distance matrix without correction 
 * (straight percent ID)
 *
 * Resolve ambiguities in sequence data by skipping
 * those nucleotides/residues
 * 
 */
void
NJ_no_correction(DMAT *dmat,
		 NJ_alignment *alignment) {

  long int i, j;
  float pdiff;

  /* compute pairwise percent identity */
  for(i=0;i<dmat->size;i++) {
    for(j=i+1;j<dmat->size;j++) {
      pdiff = 1.0 - NJ_pw_percentid(alignment, i, j);      
      dmat->val[NJ_MAP(i, j, dmat->size)] = pdiff;
    }
  }
  
  return;
}




/*
 * NJ_DNA_jc_correction() - 
 *
 * Compute the distance matrix with jukes-cantor correction
 * and assign high distance if sequence divergence exceeds
 * 0.75
 *
 *   Jukes, T.H. (1969), Evolution of protein molecules.  In H.N. Munro (Ed.),
 *   Mammalian Protein Metabolism, Volume III, Chapter 24, pp. 21-132. 
 *   New York: Academic Press
 *
 */
void
NJ_DNA_jc_correction(DMAT *dmat,
		     NJ_alignment *alignment) {
  
  long int i, j;
  long int k;
  float d, cutoff, dist;
  long int residues;

  cutoff = 0.75;
  
  for(i=0;i<dmat->size;i++) {
    for(j=i+1;j<dmat->size;j++) {
      
      k = NJ_pw_differences(alignment, i, j, &residues);
      d = 1.0 - NJ_pw_percentid(alignment, i, j);      
      
      if(d > cutoff) {
	dist = NJ_BIGDIST;
      } else {
	dist = (-0.75) * log(1.0 - (4.0/3.0)*d);
      }
      
      if(fabs(dist) < FLT_EPSILON) {
	dmat->val[NJ_MAP(i, j, dmat->size)] = 0.0;
      } else {
	dmat->val[NJ_MAP(i, j, dmat->size)] = dist;
      }
    }
  }


  
  return;
}






/*
 * NJ_PROTEIN_jc_correction() - 
 *
 * This function performs modified jukes/cantor correction on
 * a protein alignment 
 *
 *   Jukes, T.H. (1969), Evolution of protein molecules.  In H.N. Munro (Ed.),
 *   Mammalian Protein Metabolism, Volume III, Chapter 24, pp. 21-132. 
 *   New York: Academic Press
 *
 */
void
NJ_PROTEIN_jc_correction(DMAT *dmat,
			 NJ_alignment *alignment) {
  
  long int i, j;
  long int residues;
  long int diff;
  float dist, x;
  

  for(i=0;i<dmat->size;i++) {
    for(j=i+1;j<dmat->size;j++) {

      diff = NJ_pw_differences(alignment, i, j, &residues);
      
      if(!diff || !residues) {
	dist = 0.0;
      } else {

	dist = (float)diff/(float)residues;
	x = ((20.0/19.0)*dist);

	if(NJ_FLT_GT(x, 1.0)) {
	  dist = NJ_BIGDIST;
	} else {
	  dist = -(19.0/20.0) * log(1.0 - x);
	}
      }
      
      dmat->val[NJ_MAP(i, j, dmat->size)] = dist;
    }
  }
  
  return;
}






/*
 * NJ_DNA_k2p_correction() -  
 *
 * Correct a distance matrix using k2p correction using
 * cutoffs to avoid problems with logarithms.
 *
 * dist = -0.5ln(1-2P-Q) - 0.25ln(1-2Q)
 *
 * But due to the logarithms, this is only valid when
 *
 * (2P+Q <= 1)  &&  
 * (2Q <= 1)
 *
 * So assign arbitary distances when these constraints are
 * not strictly followed.
 *
 *   Kimura, M. (1980), A simple method for estimating evolutionary
 *   rates of base substitutions through comparative studies of
 *   nucleotide sequences.  J. Mol. Evol., 16, 111-120
 *
 */
void
NJ_DNA_k2p_correction(DMAT *dmat,
		      NJ_alignment *alignment) {

  long int i, j;
  float P;  /* proportion of transitions   */
  float Q;  /* proportion of transversions */
  long int nucleotides;
  long int transitions, transversions;
  float dist;
  float log_x = 0.0;  /* the params for the first log  */
  float log_y = 0.0;  /* the params for the second log */

  int blowup;   /* a flag to specify if we have a log blowup */

  
  for(i=0;i<dmat->size;i++) {
    for(j=i+1;j<dmat->size;j++) {

      blowup = 0;

      /* count the number of transitions and transversions */
      NJ_DNA_count_tt(alignment, i, j, &transitions, &transversions, &nucleotides);

      if(!nucleotides) {   /* sequences have no non-ambiguous overlap in alignment */
	P = 0.0;
	Q = 0.0;
      } else {
	P = (float)transitions   / (float)nucleotides;
	Q = (float)transversions / (float)nucleotides;
      }

      /* the first log blows up if 2*P+Q = 1.0 */
      if(NJ_FLT_EQ((2.0 * P + Q), 1.0)) {
	blowup = 1;
      } else {
	if( NJ_FLT_LT(1.0 - 2.0*P - Q, 0.0) ) {
	  blowup = 1;
	} else {
	  log_x = log(1.0 - 2.0*P - Q);
	}
      }

      /* the second log blows up if 2*Q >= 1.0 */
      if( NJ_FLT_EQ((2.0 * Q), 1.0) ||
	  NJ_FLT_GT((2.0 * Q), 1.0) ) {
	blowup = 1;
      } else {
	log_y = log(1.0 - 2.0*Q);
      }
      
      /* if our logarithms blow up, we just set the distance to the max */
      if(blowup) {
	dist = NJ_BIGDIST;
      } else {
	dist = (-0.5)*log_x - 0.25*log_y;
      }
      
      if(fabs(dist) < FLT_EPSILON) {
	dmat->val[NJ_MAP(i, j, dmat->size)] = 0.0;
      } else {
	dmat->val[NJ_MAP(i, j, dmat->size)] = dist;
      }
    }
  }
  
  return;
}




/*
 * NJ_PROTEIN_kimura_correction() - 
 *
 * Perform Kimura correction for distances derived from protein
 * alignments.
 *
 *   Kimura, M. (1983), The Neutral Theory of Molecular Evolution.
 *   p. 75., Cambridge University Press, Cambridge, England
 *
 */
void
NJ_PROTEIN_kimura_correction(DMAT *dmat,
			     NJ_alignment *alignment) {

  long int i, j;
  long int residues;
  long int diff;
  float dist;
  

  printf("NJ_PROTEIN_kimura_correction()\n");

  for(i=0;i<dmat->size;i++) {
    for(j=i+1;j<dmat->size;j++) {
      diff = NJ_pw_differences(alignment, i, j, &residues);
      
      if(!diff || !residues) {
	dist = 0.0;
      } else {
	dist = (float)diff/(float)residues;
      }
      
      if(NJ_FLT_LT(dist, 0.75)) {
	if(NJ_FLT_GT(dist, 0.0) ) {
	  dist = -log(1.0 - dist - (dist * dist/5.0) );
	}
      } else {
	if(NJ_FLT_GT(dist, 0.93) ) {
	  dist = 10.0; 
	} else {
	  dist = (float)NJ_dayhoff[ (int)((dist*1000.0)-750.0) ] / 100.0 ;
	}
      }
      
      dmat->val[NJ_MAP(i, j, dmat->size)] = dist;
    }
  }
  
  return;
}





/*
 * NJ_DNA_count_tt() - 
 *
 * Count the number of transitions and transversions
 * between two aligned DNA sequences
 *
 * This routine automatically skips ambiguities when
 * counting transitions and transversions.
 *
 */
void
NJ_DNA_count_tt(NJ_alignment *alignment,
		long int x,
		long int y,
		long int *transitions,
		long int *transversions,
		long int *residues) {

  long int tmp_transitions   = 0;
  long int tmp_transversions = 0;
  long int tmp_residues      = 0;
  char a, b;
  long int i;

  for(i=0;i<alignment->length;i++) {

    a = toupper(alignment->data[x*alignment->length+i]);
    b = toupper(alignment->data[y*alignment->length+i]);
    
    if( (a == 'A' && b == 'T') ||
	(a == 'T' && b == 'A') ||
	(a == 'A' && b == 'C') ||
	(a == 'C' && b == 'A') ||
	(a == 'T' && b == 'G') ||
	(a == 'G' && b == 'T') ||
	(a == 'C' && b == 'G') ||
	(a == 'G' && b == 'C') ) {
      tmp_transversions++;
    }
	
    if( (a == 'C' && b == 'T') ||
	(a == 'T' && b == 'C') ||
	(a == 'G' && b == 'A') ||
	(a == 'A' && b == 'G') ) {
      tmp_transitions++;
    }

    /* count the number of residues */
    if(a != NJ_AMBIGUITY_CHAR &&
       b != NJ_AMBIGUITY_CHAR ) {
      tmp_residues++;
    }

  }
  
  *transitions   = tmp_transitions;
  *transversions = tmp_transversions;
  
  if(residues) {
    *residues = tmp_residues;
  }
  
  return;
}





/*
 * NJ_pw_percentid() - 
 *
 * Given an alignment and a specification
 * for two rows, compute the pairwise
 * percent identity between the two
 *
 */
float
NJ_pw_percentid(NJ_alignment *alignment,
		long int x,
		long int y) {
  
  float pid;
  long int i;
  long int residues;
  long int same;
  char c1, c2;

  residues = 0;
  same     = 0;
  for(i=0;i<alignment->length;i++) {

    c1 = alignment->data[x*alignment->length+i];
    c2 = alignment->data[y*alignment->length+i];
    
    if( c1 != NJ_AMBIGUITY_CHAR ||
        c2 != NJ_AMBIGUITY_CHAR ) {
      
      residues++;

      if(c1 == c2) {
	same++;
      }
    }

  }

  pid = (float)same/(float)residues;
  
  return(pid);
}



/*
 * NJ_pw_differences() - 
 *
 * Given an alignment and a specification
 * for two rows in the alignment, compute the
 * number of differences between the two sequences
 *
 * With respect to ambiguity codes, we will want to 
 * disregard those sites entirely in our count.
 *
 */
long int
NJ_pw_differences(NJ_alignment *alignment,
		  long int x,
		  long int y,
		  long int *residues) {

  long int i;
  long int diff;
  char c1, c2;
  long int tmp_residues;
  
  diff         = 0;
  tmp_residues = 0;
  for(i=0;i<alignment->length;i++) {

    c1 = alignment->data[x*alignment->length+i];
    c2 = alignment->data[y*alignment->length+i];
    
    if( c1 != NJ_AMBIGUITY_CHAR ||
        c2 != NJ_AMBIGUITY_CHAR ) {
      
      tmp_residues++;

      if(c1 != c2) {
	diff++;
      }
    }

  }

  *residues = tmp_residues;

  return(diff);
}







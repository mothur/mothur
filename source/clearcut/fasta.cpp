/*
 * fasta.c
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
 * Functions for parsing FASTA formatted alignment files
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
#include <string.h>
#include <ctype.h>

#include "clearcut.h"
#include "common.h"
#include "fasta.h"


#define NJ_NUM_DNA_AMBIGUITY_SYMS 14
static const char NJ_dna_ambiguity_syms[NJ_NUM_DNA_AMBIGUITY_SYMS] = 
{
  'M', 'R', 'W', 'S', 'Y', 'K',
  'V', 'H', 'D', 'B', 'X', 'N',
  '-', '.'
};


#define NJ_NUM_PROTEIN_AMBIGUITY_SYMS 6
static const char NJ_protein_ambiguity_syms[NJ_NUM_PROTEIN_AMBIGUITY_SYMS] =
{
  'X', 'B', 'Z', '*', '-', '.'
};

#define NJ_NUM_DNA_SYMS 5
static const char NJ_dna_syms[NJ_NUM_DNA_SYMS] = 
{
  'A', 'G', 'C', 'T', 'U'
};


#define NJ_NUM_PROTEIN_SYMS 20
static const char NJ_protein_syms[NJ_NUM_PROTEIN_SYMS] = 
{
  'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
};





/*
 * NJ_is_whitespace() - Check to see if character is whitespace
 *
 * INPUTS:
 * -------
 *     c -- character to check
 *
 * RETURNS:
 * --------
 *    int -- 0 if not whitespace
 *           1 if whitespace
 */
static inline
int
NJ_is_whitespace(char c) {
  if( c == ' '  ||   /* space           */
      c == '\n' ||   /* newline         */
      c == '\r' ||   /* carriage-return */
      c == '\v' ||   /* vertical tab    */
      c == '\f' ||   /* form feed       */
      c == '\t' ) {  /* horizontal tab  */
    return(1);
  } else {
    return(0);
  }
}





/*
 * NJ_is_dna() - 
 *
 * Determines if the given symbol is DNA
 *
 * RETURNS: 1 if DNA
 *          0 if not DNA
 *
 */
static inline
int
NJ_is_dna(char c) {

  int i;
  char up_c;
  
  up_c = toupper(c);
  
  for(i=0;i<NJ_NUM_DNA_SYMS;i++) {
    if(up_c == NJ_dna_syms[i]) {
      return(1);
    }
  }

  return(0);
}





/*
 * NJ_is_dna_ambiguity() - 
 *
 * Determines if the given symbol is a 
 * DNA ambiguity code
 *
 * RETURNS: 1 if DNA Ambiguity Code
 *          0 if not DNA Ambiguity Code
 *
 */
static inline
int
NJ_is_dna_ambiguity(char c) {
  
  int i;
  char up_c;
  
  up_c = toupper(c);
  
  for(i=0;i<NJ_NUM_DNA_AMBIGUITY_SYMS;i++) {
    if(up_c == NJ_dna_ambiguity_syms[i]) {
      return(1);
    }
  }
  
  return(0);
}



/*
 * NJ_is_protein() - 
 *
 * Determines if supplied symbol is amino acid symbol
 *
 * RETURNS: 1 if amino acid
 *          0 if not amino acid
 *
 */
static inline
int
NJ_is_protein(char c) {

  int i;
  char up_c;
  
  up_c = toupper(c);
  
  for(i=0;i<NJ_NUM_PROTEIN_SYMS;i++) {
    if(up_c == NJ_protein_syms[i]) {
      return(1);
    }
  }
  
  return(0);
}




/*
 * NJ_is_protein_ambiguity() - 
 *
 * Determines if supplied symbol is amino acid ambiguity symbol
 *
 * RETURNS: 1 if amino acid ambiguity symbol
 *          0 if not amino acid ambiguity symbol
 *
 */
static inline
int 
NJ_is_protein_ambiguity(char c) {

  int i;
  char up_c;
  
  up_c = toupper(c);

  for(i=0;i<NJ_NUM_PROTEIN_AMBIGUITY_SYMS;i++) {
    if(up_c == NJ_protein_ambiguity_syms[i]) {
      return(1);
    }
  }
  
  return(0);
}






/*
 * NJ_read_fasta() - A function for inputing FASTA sequences into an alignment
 *                   data structure
 *
 *
 * INPUTS:
 * -------
 *   nj_args -- A pointer to a data structure containing command-line args
 *
 * RETURNS:
 * --------
 *   NJ_alignment * -- A pointer to a multiple sequence alignment
 *
 * DESCRIPTION:
 * ------------
 *
 * 
 * This function implements a state machine parser for parsing FASTA-formatted
 * multiple sequence alignment files.  
 * 
 * Example Input:
 *
 *   > sequence1
 *   ATAGATATAGATTAGAATAT----TATAGATAT----ATATAT-TTT-
 *   > sequence2 
 *   --ATAGATA---ATATATATATTTT--GTCTCATAGT---ATATGCTT
 *   > sequence3
 *   TTATAGATA---ATATATATATTTTAAGTCTCATAGT-A-ATATGC--
 * 
 * This function will parse alignments for DNA or protein, and will do
 * so mindful of ambiguity codes for these kinds of sequences.  All 
 * ambiguity codes are ignored by this program for the purposes of 
 * computing a distance matrix from a multiple alignment.  By design, 
 * this program does not auto-detect DNA vs. Protein, and requires that 
 * the user explictly specify that on the command-line.
 *
 * Gaps can be represented either with the '-' or '.' characters.
 * 
 * Newlines and other whitespace are allowed to be interspersed 
 * throughout the sequences.
 *
 * Taxon labels are required to be unique, and they must start with 
 * an alphabetic character (not a number, etc.).  The parser will read
 * the first token after the > character in the description line up until the
 * first whitespace and use that for the taxon label.
 *
 * For example, in the line "> taxon1 is homo sapien", the taxon label will be 
 * "taxon1"
 *
 */
NJ_alignment *
NJ_read_fasta(NJ_ARGS *nj_args) {

  FILE *fp  = nullptr;
  char *buf = nullptr;
  char *ptr = nullptr;
  NJ_alignment *alignment = nullptr;

  char c;
  int state;
  long int index, x, seq;
  long int i;
  long int bufsize, nseqs = NJ_INITIAL_NSEQS;
  int first_sequence_flag;


  

  /* 
   * In this function, we implement a FASTA alignment parser which
   * reads in an alignment character-by-character, maintaining state
   * information which guides the parser.
   *
   * The program reads either DNA or Protein alignments.  All title lines
   * and sequences can be arbitrarily long.  Gaps can be represented by 
   * "-" or "." characters.  
   *
   * Ambiguity codes are also handled.
   * 
   */

  /* 
   * We can't handle reading fasta input unless the user explicity 
   * specifies the input type...just to be sure.
   */
  if( (!nj_args->dna_flag && !nj_args->protein_flag) ||
      (nj_args->dna_flag  &&  nj_args->protein_flag) ) {
    fprintf(stderr, "Clearcut: Explicitly specify protein or DNA\n");
    goto XIT_BAD;
  }

  /* open specified fasta file here */
  if(nj_args->stdin_flag) {
    fp = stdin;
  } else {
    fp = fopen(nj_args->infilename, "r");
    if(!fp) {
      fprintf(stderr, "Clearcut: Failed to open input FASTA file: %s\n", nj_args->infilename);
      perror("Clearcut");
      goto XIT_BAD;
    }
  }

  /* allocate the initial buffer */
  bufsize = NJ_INITIAL_BUFFER_SIZE;
  buf = (char *)calloc(bufsize, sizeof(char));
  
  /* allocate the alignment container here */
  alignment = (NJ_alignment *)calloc(1, sizeof(NJ_alignment));
  
  /* allocate initial title array */
  //  printf("allocating initial title array\n");
  alignment->titles = (char **)calloc(NJ_INITIAL_NSEQS, sizeof(char *));

  /* make sure that we successfully allocated memory */
  if(!buf || !alignment || !alignment->titles) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_read_fasta()\n");
    goto XIT_BAD;
  }

  /* a flag */
  first_sequence_flag = 1;  
  
  index  = 0;  /* tracks the position in buffer     */
  x      = 0;  /* tracks the position on sequence   */
  seq    = 0;  /* tracks the active sequence        */

  /* intitial state of state machine */
  state  = NJ_FASTA_MODE_UNKNOWN;

  while(1) {
    
    /* get the next character */
    c = fgetc(fp);
    if(feof(fp)) {

      if(state == NJ_FASTA_MODE_SEQUENCE) {
	buf[index+1] = '\0';

	/* copy buf to alignment */
	for(i=1;i<=alignment->length;i++) {
	  alignment->data[seq*alignment->length+i-1] = buf[i];
	}
      }

      break;
    }

    /* make sure our dynamic buffer is big enough */
    if(index >= bufsize) {
      bufsize *= 2;
      buf = (char *)realloc(buf, bufsize);
      if(!buf) {
	fprintf(stderr, "Clearcut: Memory allocation error in NJ_read_fasta()\n");
	goto XIT_BAD;
      }
    }

    switch(state) {
      
    case NJ_FASTA_MODE_UNKNOWN:
      
      if(!NJ_is_whitespace(c)) {
	if(c == '>') {
	  state = NJ_FASTA_MODE_TITLE;
	  buf[0] = '>';
	} else {
	  goto XIT_BAD;
	}
      }

      break;

    case NJ_FASTA_MODE_TITLE:

      if( c == '\n' ||
	  c == '\r' ) {

	buf[index] = '\0';
	state = NJ_FASTA_MODE_SEQUENCE;
	index = 0;
	x     = -1;

	/* make sure we've allocated enough space for titles and sequences */
	if(seq == nseqs) {

	  //	  printf("realloc().  seq = %d, nseqs = %d\n", seq, nseqs);

	  nseqs *= 2;

	  alignment->titles = (char **)realloc(alignment->titles, nseqs*sizeof(char *));
	  if(!alignment->titles) {
	    fprintf(stderr, "Clearcut:  Memory allocation error in NJ_read_fasta()\n");
	    goto XIT_BAD;
	  }

	  alignment->data = (char *)realloc(alignment->data, alignment->length*nseqs*sizeof(char));
	  if(!alignment->data) {
	    fprintf(stderr, "Clearcut: Allocation error in NJ_read_fasta()\n");
	    goto XIT_BAD;
	  }
	}

	// printf("Allocating %d bytes for title %d: %s\n", (int)strlen(buf), (int)seq, buf);

	alignment->titles[seq] = (char *)calloc(strlen(buf), sizeof(char));
	if(!alignment->titles[seq]) {
	  fprintf(stderr, "Clearcut:  Memory allocation error in NJ_read_fasta()\n");
	  goto XIT_BAD;
	}

	/* lets forward to the first non-space (space/tab) character after the '>' */

	if(first_sequence_flag) {
	  ptr = buf;
	} else {
	  ptr = &buf[1];
	}
	while(*ptr == '\t' || *ptr == ' ') {
	  ptr++;
	}
	sscanf(ptr, "%s", alignment->titles[seq]);  /* get the first word and use as the title */

	alignment->nseq++;
      }
	
      buf[index++] = c;

      break;


    case NJ_FASTA_MODE_SEQUENCE:

      if(c == '>') {

	if(first_sequence_flag) {
	  first_sequence_flag = 0;

	  /* allocate our alignment data section here */
	  alignment->length = index-1;

	  nseqs = NJ_INITIAL_NSEQS;
	  alignment->data = (char *)calloc(alignment->length*nseqs, sizeof(char));
	  if(!alignment->data) {
	    fprintf(stderr, "Clearcut: Allocation error in NJ_read_fasta()\n");
	    goto XIT_BAD;
	  }
	} 
	
	if(!first_sequence_flag) {
	  if(index-1 < alignment->length) {
	    fprintf(stderr, "Clearcut: Sequences must be of uniform length in alignment at sequence %ld\n", seq);
	    goto XIT_BAD;
	  }
	}

	/* re-allocate if necessary */
	/*
	if(seq >= nseqs) {
	  nseqs *= 2;
	  alignment->data = (char *)realloc(alignment->data, alignment->length*nseqs*sizeof(char));
	  if(!alignment->data) {
	    fprintf(stderr, "Clearcut: Allocation error in NJ_read_fasta()\n");
	    goto XIT_BAD;
	  }
	}
	*/

	/* copy buf to alignment */
	for(i=1;i<=alignment->length;i++) {
	  alignment->data[seq*alignment->length+i-1] = buf[i];
	}
	  
	state = NJ_FASTA_MODE_TITLE;
	index = 1;
	x     = 1;
	  
	buf[0] = c;
	
	seq++;

      } else {
	  
	if(NJ_is_whitespace(c)) {
	  break;
	}
		
	if(!first_sequence_flag) {
	  if(index-1 >= alignment->length) {
	    fprintf(stderr, "Clearcut: Sequences must be of uniform length in alignment at sequence %ld\n", seq);
	    goto XIT_BAD;
	  }
	}


	/* 
	 * Here we check to make sure that the symbol read is appropriate
	 * for the type of data the user specified.  (dna or protein).
	 * We also handle ambiguity codes by converting them to a specific
	 * assigned ambiguity code character.  Ambiguity codes are ignored
	 * when computing distances
	 */

	if(nj_args->dna_flag) {
	  if(NJ_is_dna(c)) {
	    buf[index++] = toupper(c);
	  } else {
	    if(NJ_is_dna_ambiguity(c)) {
	      buf[index++] = NJ_AMBIGUITY_CHAR;
	    } else {
	      fprintf(stderr, "Clearcut: Unknown symbol '%c' in nucleotide sequence %ld.\n", c, seq);
	      goto XIT_BAD;
	    }
	  }
	} else if(nj_args->protein_flag) {
	  if(NJ_is_protein(c)) {
	    buf[index++] = toupper(c);
	  } else {
	    if(NJ_is_protein_ambiguity(c)) {
	      buf[index++] = NJ_AMBIGUITY_CHAR;
	    } else {
	      fprintf(stderr, "Clearcut: Unknown symbol '%c' in protein sequence %ld.\n", c, seq);
	      goto XIT_BAD;
	    }
	  }
	}

      }

      break;
	
    default:
      goto XIT_BAD;
	
      break;

    }
  }
  
  if(index-1 != alignment->length) {
    fprintf(stderr, "Clearcut: Sequences must be of uniform length in alignment at sequence %ld\n", seq);
    goto XIT_BAD;
  }
  
  /* check for duplicate taxon labels */
  if(!NJ_taxaname_unique(alignment)) {
    goto XIT_BAD;
  }
  
  return(alignment);

  
 XIT_BAD:

  if(fp) {
    fprintf(stderr, "Clearcut: Fatal error parsing FASTA file at file offset %ld.\n", ftell(fp));
  }
  
  if(buf) {
    free(buf);
  }
  
  NJ_free_alignment(alignment);

  return(nullptr);
}




/*
 * NJ_print_alignment() - Print multiple sequence alignment (for debugging)
 *
 * INPUTS:
 * -------
 *    alignment -- A pointer to the alignment
 *
 * RETURNS:
 * --------
 *    NONE
 *
 */
void
NJ_print_alignment(NJ_alignment *alignment) {
  
  long int i, j;
  
  printf("nseq = %ld, length = %ld\n", alignment->nseq, alignment->length);
  
  for(i=0;i<alignment->nseq;i++) {
    
    printf("> %s\n", alignment->titles[i]);

    for(j=0;j<alignment->length;j++) {
      printf("%c", alignment->data[i*alignment->length+j]);
    }

    printf("\n");
  }

  return;
}







/*
 *
 * NJ_free_alignment() - Free all of the memory allocated for the
 *                       multiple sequence alignment 
 *
 * INPUTS:
 * -------
 *   alignment -- A pointer to the multiple sequence alignment
 *
 * RETURNS:
 * --------
 *    NONE
 *
 */
void
NJ_free_alignment(NJ_alignment *alignment) {
  
  long int i;
  
  if(alignment) {

    /* free the allocated titles */
    if(alignment->titles) {
      for(i=0;i<alignment->nseq;i++) {

	if(alignment->titles[i]) {
	  free(alignment->titles[i]);
	}
      }
      
      free(alignment->titles);
    }

    /* free the alignment data */
    if(alignment->data) {
      free(alignment->data);
    }

    /* free the alignment itself */
    free(alignment);
  }

  return;
}




/*
 * NJ_taxaname_unique() - Check to see if taxanames are unique in alignment
 *
 * INPUTS:
 * -------
 *  alignment -- a pointer to a multiple sequence alignment
 *
 * OUTPUTS:
 * --------
 *  int -- 0 if all taxanames in alignment are unique
 *         1 if all taxanames in alignment are NOT unique
 *
 *
 * DESCRIPTION:
 * ------------
 *
 * Check to see if the taxanames in the alignment are unique.  It
 * will be impossible to make sense of the final tree if the taxon
 * labels are not unqiue.
 *
 */
int
NJ_taxaname_unique(NJ_alignment *alignment) {
  
  long int i, j;
  
  for(i=0;i<alignment->nseq;i++) {
    for(j=i+1;j<alignment->nseq;j++) {
        if(!strcmp(alignment->titles[i], 
		 alignment->titles[j])) {
	fprintf(stderr, "Clearcut: Taxa %ld and %ld (%s) do not have unique taxon labels.\n", 
		i, j, alignment->titles[i]);
	return(0);
      }
    }
  }
  
  return(1);
}


void
NJ_print_titles(NJ_alignment *alignment) {

  int i;

  for(i=0;i<alignment->nseq;i++) {
    printf("%d: %s\n", i, alignment->titles[i]);
  }

  return;
}

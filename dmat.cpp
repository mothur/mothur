/*
 * dmat.c
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
 * Distance matrix parser
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
#include <float.h>
#include <errno.h>
#include <string.h>



#include "common.h"
#include "clearcut.h"

#include "dmat.h"





/*
 *
 * NJ_is_alpha() - determine if character is an alphabetic character
 *
 * INPUT:
 * ------
 *  c -- character to test
 *
 * RETURN:
 * -------
 *   int -- 1 if character is alphabetic (A-Z || a-z)
 *          0 if character is NOT alphabetic
 *
 */
static inline
int 
NJ_is_alpha(char c) {

  if( (c >= 'A' && c <= 'Z') ||
      (c >= 'a' && c <= 'z') ) {
    return(1);
  } else {
    return(0);
  }
}



/*
 *
 * NJ_is_whitespace() - determine if character is a whitespace character
 *
 * INPUT:
 * ------
 *  c -- character to test
 *
 * RETURN:
 * -------
 *   int -- 1 if character is whitespace (space, tab, CR, LF)
 *          0 if character is NOT whitespace
 *
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
 *
 * NJ_is_number() - determine if character is a number
 *
 * INPUT:
 * ------
 *  c -- character to test
 *
 * RETURN:
 * -------
 *   int -- 1 if character is a number (0-9)
 *          0 if character is NOT a number
 *
 */
static inline
int
NJ_is_number(char c) {
  if(c >= '0' && c <= '9') {
    return(1);
  } else {
    return(0);
  }
}



/*
 * NJ_is_distance() - check if string is a properly formatted distance value
 *
 */
static inline
int
NJ_is_distance(char *token) {
  
  int i;
  char c;
  int exponent_state;
  int expsign_state;
  int dpoint_state;

  /* if token is NULL return failure */
  if(!token) {
    return(0);
  }
  
  exponent_state = 0;
  expsign_state  = 0;
  dpoint_state   = 0;
  
  /* The first character must be a number, a decimal point or a sign */
  c = token[0];
  if(!NJ_is_number(c) &&
     c != '.'         &&
     c != '-'         &&
     c != '+' )  {
    
    goto BAD;
  } 
  
  /* 
   * if the first character is not a number, and string is only one 
   * character long, then we return failure.
   */
  if(strlen(token) == 1) {
    if(!NJ_is_number(c)) {
      goto BAD;
    }
  }
  
  for(i=0;i<strlen(token);i++) {

    c = token[i];
    
    /* make sure that all chars in dist string are in list of valid chars */
    if(!NJ_is_number(c) &&
       c != '.'         &&
       c != '-'         &&
       c != '+'         &&
       c != 'e'         &&
       c != 'E'         ) {

      goto BAD;
    }

    /* not the first char and we are not in exponent state but read (+,-) */
    if(i>0 && !exponent_state) {
      if(c == '-' || 
	 c == '+') {
	goto BAD;
      }
    }

    /* if we are in the exponent state, and we've already seen a sign */
    if(exponent_state && expsign_state) {
      if(c == '-' ||
	 c == '+') {
	goto BAD;
      }
    }
    
    /* if we are in the exponent state and we see a decimal point */
    if(exponent_state) {
      if(c == '.') {
	goto BAD;
      }
    }
    
    /* if we are in the exponent state and see another e or E */
    if(exponent_state) {
      if(c == 'e' ||
	 c == 'E') {
	goto BAD;
      }
    }
    
    /* if we are dpoint_state and see another decimal point */
    if(dpoint_state) {
      if(c == '.') {
	goto BAD;
      }
    }
    
    
    /* enter the exponent state if we need to */
    if(!exponent_state) {
      if(c == 'e' ||
	 c == 'E') {
	exponent_state = 1;
      }
    }

    /* enter the expsign_state if we need to */
    if(exponent_state && !expsign_state) {
      if(c == '-' ||
	 c == '+') {
	expsign_state = 1;
      }
    }

    /* if not in dpoint state and we see a dpoint */
    if(!dpoint_state) {
      if(c == '.') {
	dpoint_state = 1;
      }
    }

  }
  
  /* the token must end in a number char */
  if(!NJ_is_number(token[strlen(token)-1])) {
    goto BAD;
  }
  
  /* token is a valid numerical distance */
  return(1);
  
 BAD:

  /* token is invalid distance format */
  return(0);
}




/*
 * NJ_is_label() - 
 *
 * Simply, if token is not a valid number, then it is a name
 *
 */
static inline
int
NJ_is_label(char *token) {
  if(NJ_is_distance(token)) {
    return(0);
  } else {
    return(1);
  }
}



/*
 * NJ_get_token() - get a token from an input stream 
 *
 */
static inline
int
NJ_get_token(FILE *fp,
	     NJ_DIST_TOKEN *token) {

  char c;
  int index;

  c = fgetc(fp);
  if(feof(fp)) {
    token->type = NJ_EOF_STATE;
    return(token->type);
  }

  if(NJ_is_whitespace(c)) {
    token->buf[0] = c;
    token->buf[1] = '\0';
    token->type   = NJ_WS_STATE;

    return NJ_WS_STATE;
  } 

  index = 0;
  while(!NJ_is_whitespace(c)) {

    /* reallocate our buffer if necessary */
    if(index >= token->bufsize) {
      token->bufsize *= 2;
      token->buf = (char *)realloc(token->buf, token->bufsize*sizeof(char));
      if(!token->buf) {
	fprintf(stderr, "Clearcut: Memory allocation error in NJ_get_token()\n");
	exit(-1);
      }
    }

    token->buf[index++] = c;
    
    c = fgetc(fp);
    if(feof(fp)) {
      token->type = NJ_EOF_STATE;
      break;
    }
  }
  
  token->buf[index] = '\0';
  
  if(token->type != NJ_EOF_STATE) {

    if(NJ_is_distance(token->buf)) {
      token->type = NJ_FLOAT_STATE;
    } else {
      token->type = NJ_NAME_STATE;
    }

  }
  
  return(token->type);
}




 

/* 
 * NJ_parse_distance_matrix() -- Takes a filename and returns a distance matrix
 *
 *
 * INPUT:
 * ------
 *   nj_args -- a pointer to a structure containing the command-line arguments
 *
 * OUTPUT:
 * -------
 *   <DMAT *> -- NULL  (failure)
 *            -- A pointer to a populated distance matrix
 *
 * DESCRIPTION:
 * ------------
 *   This function implements a simple state machine to parse a distance matrix
 *   in approximate PHYLIP format.  This function auto-detects whether the 
 *   distance matrix is in upper, lower, or fully-symmetric format and handles
 *   it accordingly.  For full/symmetric matrices, values must be symmetric 
 *   around the diagonal, which is required to be zeroes.  Names and values must
 *   be separated by whitespace (space, tab, newlines, etc.).  Taxon labels can
 *   include numbers, but must start with non-numerical symbols.
 *
 *
 *   *** UPPER FORMAT EXAMPLE ***
 * 
 *   4
 *   seq1 0.2 0.3 0.1
 *   seq2     0.2 0.3
 *   seq3         0.1
 *   seq4
 * 
 *   *** LOWER FORMAT EXAMPLE ***
 *
 *   4
 *   seq1
 *   seq2 0.3
 *   seq3 0.2 0.4
 *   seq4 0.3 0.1 0.3
 *
 *   *** SYMMETRIC (FULL) EXAMPLE ***
 *  
 *   4
 *   seq1 0.0 0.3 0.5 0.3
 *   seq2 0.3 0.0 0.1 0.2
 *   seq3 0.5 0.1 0.0 0.9
 *   seq4 0.3 0.2 0.9 0.0
 *
 *  Values in the distance matrix can be positive or negative, integers or
 *  real values.  Values can also be parsed in exponential notation form.
 * 
 */
DMAT *
NJ_parse_distance_matrix(NJ_ARGS *nj_args) {
  
  DMAT *dmat           = NULL;
  FILE *fp            = NULL;
  NJ_DIST_TOKEN *token = NULL;

  int state, dmat_type;
  int row;
  int fltcnt;
  int x, y, i;
  int numvalread;
  int expectedvalues = -1;
  float val;
  int first_state = 0;


  /* allocate our distance matrix and token structure */
  dmat = (DMAT *)calloc(1, sizeof(DMAT));
  token = (NJ_DIST_TOKEN *)calloc(1, sizeof(NJ_DIST_TOKEN));
  if(token) {
    token->bufsize = NJ_INITIAL_BUFSIZE;
    token->buf     = (char *)calloc(token->bufsize, sizeof(char));
  }
  if(!dmat || !token || !token->buf) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_parse_distance_matrix()\n");
    goto XIT_BAD;
  }
  
  /* open distance matrix file here */
  if(nj_args->stdin_flag) {
    fp = stdin;
  } else {
    fp = fopen(nj_args->infilename, "r");
    if(fp==NULL) {
      fprintf(stderr, "Clearcut: Could not open distance matrix: %s\n", nj_args->infilename);
      perror("Clearcut");
      goto XIT_BAD;
    }
  }

  /* get the number of taxa in this file */
  fscanf(fp, "%ld\n", &dmat->ntaxa);
  if(dmat->ntaxa < 2) {
    fprintf(stderr, "Clearcut: Invalid number of taxa in distance matrix\n");

    goto XIT_BAD;
  }

  /* set our initial working size according to the # of taxa */
  dmat->size = dmat->ntaxa;

  /* allocate space for the distance matrix values here */
  dmat->val = 
    (float *)calloc(NJ_NCELLS(dmat->ntaxa), sizeof(float));
  if(!dmat->val) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_parse_distance_matrix()\n");
    goto XIT_BAD;
  }

  /*  taxa names */
  dmat->taxaname = (char **)calloc(dmat->ntaxa, sizeof(char *));
  if(!dmat->taxaname) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_parse_distance_matrix()\n");
    goto XIT_BAD;
  }

  /* set the initial state of our state machine */
  dmat_type   = NJ_PARSE_UNKNOWN;
  row         = -1;
  fltcnt      = 0;
  numvalread  = 0;


  /* read the input one character at a time to drive simple state machine */
  state = NJ_get_token(fp, token);
  while(state != NJ_EOF_STATE) {
    
    switch(state) {

    case NJ_NAME_STATE:
      
      if(first_state == 0) {
	first_state = 1;
      }

      row++;

      if(row > 0 && dmat_type == NJ_PARSE_UNKNOWN) {

	if(fltcnt == dmat->ntaxa) {

	  dmat_type = NJ_PARSE_SYMMETRIC;
	  expectedvalues = dmat->ntaxa * dmat->ntaxa;

	} else if (fltcnt == dmat->ntaxa-1) {

	  dmat_type = NJ_PARSE_UPPER;
	  expectedvalues = ((dmat->ntaxa) * (dmat->ntaxa-1)) / 2;

	  /* shift everything in first row by one char */
	  for(i=dmat->ntaxa-2;i>=0;i--) {
	    dmat->val[i+1] = dmat->val[i];
	  }

	} else if (fltcnt == 0) {

	  dmat_type = NJ_PARSE_LOWER;
	  expectedvalues = ((dmat->ntaxa) * (dmat->ntaxa-1)) / 2;

	} else {
	  goto XIT_BAD;
	}
      }
      
      if(row >= dmat->ntaxa) {
	goto XIT_BAD;
      }
      
      /* allocate space for this taxon label */
      dmat->taxaname[row] = (char *)calloc(strlen(token->buf)+1, sizeof(char));
      if(!dmat->taxaname[row]) {
	fprintf(stderr, "Clearcut: Memory allocation error in NJ_parse_distance_matrix()\n");
	goto XIT_BAD;
      }
      
      strcpy(dmat->taxaname[row], token->buf); 
      
      fltcnt = 0;

      break;


    case NJ_FLOAT_STATE:

      if(first_state == 0) {
	goto XIT_BAD;
      }

      val = atof(token->buf);
      if(errno) {
	fprintf(stderr, "Clearcut: Distance value out-of-range.\n");
	goto XIT_BAD;
      }
      
      x = row;
      y = fltcnt;

      switch(dmat_type) {

      case NJ_PARSE_UNKNOWN:

	dmat->val[NJ_MAP(x, y, dmat->size)] = val;

	break;

      case NJ_PARSE_SYMMETRIC:
	
	if(fltcnt >= dmat->ntaxa) {
	  fprintf(stderr, "Clearcut: Incorrect number of distance values on row.\n");
	  goto XIT_BAD;
	}

	if(x < y) {
	  dmat->val[NJ_MAP(x, y, dmat->size)] = val;
	} else if(x > y) {
	  if(!NJ_FLT_EQ(val, dmat->val[NJ_MAP(y, x, dmat->size)])) {
	    fprintf(stderr, "Clearcut: Full matrices must be symmetric.\n");
	    goto XIT_BAD;
	  }
	} else {
	  if(!NJ_FLT_EQ(val, 0.0)) {
	    fprintf(stderr, "Clearcut: Values along the diagonal in a symmetric matrix must be zero.\n");
	    goto XIT_BAD;

            }
	}

	break;
	
      case NJ_PARSE_UPPER:

	if(fltcnt > dmat->ntaxa-row) {
	  fprintf(stderr, "Clearcut: Incorrect number of distance values on row.\n");
	  goto XIT_BAD;
	}
	
	dmat->val[NJ_MAP(x, x+y+1, dmat->size)] = val;
	
	break;

      case NJ_PARSE_LOWER:
	
	if(fltcnt > row-1) {
	  fprintf(stderr, "Clearcut: Incorrect number of distance values on row.\n");
	  goto XIT_BAD;
	}

	dmat->val[NJ_MAP(y, x, dmat->size)] = val;
	
	break;
	
      default:
	
	goto XIT_BAD;
	
	break;
      }

      fltcnt++;
      numvalread++;
      
      break;

    case NJ_WS_STATE:

      break;

    case NJ_EOF_STATE:

      if(first_state == 0) {
	goto XIT_BAD;
      }

      break;

    default:

      fprintf(stderr, "Clearcut: Unknown state in distance matrix parser.\n");
      break;

    }

    /* get next token from stream */
    state = NJ_get_token(fp, token);
  }


  /* 
   * At the end, if we have not read the number of values that we predicted
   * we would need, then there was a problem and we need to punt.
   */
  if(numvalread != expectedvalues) {
    fprintf(stderr, "Clearcut: Incorrect number of values in the distance matrix.\n");
    goto XIT_BAD;
  }
  
  /* special check to make sure first value read is 0.0 */
  if(dmat_type == NJ_PARSE_SYMMETRIC) {
    if(!NJ_FLT_EQ(dmat->val[NJ_MAP(0, 0, dmat->size)], 0.0)) {
      fprintf(stderr, "Clearcut: Values along the diagonal in a symmetric matrix must be zero.\n");
      goto XIT_BAD;
    }
  }

  
  /* now lets allocate space for the r and r2 columns */
  dmat->r  = (float *)calloc(dmat->ntaxa, sizeof(float));
  dmat->r2 = (float *)calloc(dmat->ntaxa, sizeof(float));
  if(!dmat->r || !dmat->r2) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_parse_distance_matrix()\n");
    goto XIT_BAD;
  }
  
  /* track some memory addresses */
  dmat->rhandle   = dmat->r;
  dmat->r2handle  = dmat->r2;
  dmat->valhandle = dmat->val;

  /* close matrix file here */
  if(!nj_args->stdin_flag) {
    fclose(fp);
  }
  
  if(token) {
    if(token->buf) {
      free(token->buf);
    }
    free(token);
  }

  return(dmat);



  /* clean up our partial progress */
 XIT_BAD:

  if(fp) {
    fprintf(stderr, "Clearcut: Syntax error in distance matrix at offset %ld.\n", ftell(fp));
  }

  /* close matrix file here */
  if(!nj_args->stdin_flag) {
    if(fp) {
      fclose(fp);
    }
  }

  /* if we have a valid dmat (partial or complete), we need to free it */
  if(dmat) {
    NJ_free_dmat(dmat);
  }
  
  if(token) {
    if(token->buf) {
      free(token->buf);
    }
    free(token);
  }
  
  return(NULL);
}







/*
 * NJ_output_matrix() - Output a distance matrix to the specified file 
 *
 *
 * INPUTS:
 * -------
 *  nj_args -- a pointer to a data structure holding the command-line args
 *     dmat -- a pointer to a distance matrix
 *
 *
 * RETURNS:
 * --------
 *   NOTHING
 *
 *
 * DESCRIPTION:
 * ------------
 *   If the appropriate flag was specified in the command-line, this function
 *   now outputs the parsed or computed distance matrix to a file.  This 
 *   can be useful if generating a distance matrix was the primary goal of 
 *   running the program, or if one wanted to debug and/or verify the
 *   correctness of the program.
 *
 *   Currently this function outputs full/symmetric matrices only.
 *
 */
void
NJ_output_matrix(NJ_ARGS *nj_args,
		 DMAT *dmat) {
  
  FILE *fp = NULL;
  long int i, j;


  
  /* if we haven't specieid matrixout, return immediately */
  if(!nj_args->matrixout) {
    return;
  }
  
  /* open the specified matrix file for writing */
  fp = fopen(nj_args->matrixout, "w");
  if(!fp) {
    fprintf(stderr, "Clearcut: Could not open matrix file %s for output.\n", nj_args->matrixout);
    return;
  }

  /* output the number of taxa in the matrix */
  fprintf(fp, "   %ld\n", dmat->size);

  fprintf(fp, "%s\n", dmat->taxaname[0]); // print the first taxon name outside of the main loop

  for(i=1;i<dmat->size;i++) {
    
    /* output taxaname */
    fprintf(fp, "%s\t", dmat->taxaname[i]);

    for(j=0;j<i;j++) {
      if(nj_args->expdist) {  /* exponential notation (or not) */
	fprintf(fp, "%e ", dmat->val[NJ_MAP(j,i,dmat->size)]);  
      } else {
	fprintf(fp, "%f ", dmat->val[NJ_MAP(j,i,dmat->size)]);
      }
    }
    
    fprintf(fp, "\n");
  }

#ifdef FULL_SYMMETRIC_MATRIX 

  /* output the number of taxa in the matrix */
  fprintf(fp, "   %ld\n", dmat->size);
  for(i=0;i<dmat->size;i++) {
    
    /* output taxaname */
    fprintf(fp, "%s\t", dmat->taxaname[i]);

    for(j=0;j<dmat->size;j++) {
      if(i>j) {
	if(nj_args->expdist) {  /* exponential notation (or not) */
	  fprintf(fp, "%e ", dmat->val[NJ_MAP(j,i,dmat->size)]);  
	} else {
	  fprintf(fp, "%f ", dmat->val[NJ_MAP(j,i,dmat->size)]);
	}
      } else if(i<j) {
	if(nj_args->expdist) {  /* exponential notation (or not) */
	  fprintf(fp, "%e ", dmat->val[NJ_MAP(i,j,dmat->size)]);
	} else {
	  fprintf(fp, "%f ", dmat->val[NJ_MAP(i,j,dmat->size)]);
	}
      } else {
	if(nj_args->expdist) {  /* exponential notation (or not) */
	  fprintf(fp, "%e ", 0.0);
	} else {
	  fprintf(fp, "%f ", 0.0);
	}
      }
    }
    
    fprintf(fp, "\n");
  }

#endif // FULL_SYMMETRIC_MATRIX
  
  /* close the file here */
  if(fp) {
    fclose(fp);
  }
  
  return;
}






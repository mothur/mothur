
/*
 * clearcut.c
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
 * An implementation of the Relaxed Neighbor-Joining algorithm 
 *  of Evans, J., Sheneman, L., and Foster, J.
 *
 *
 * AUTHOR:
 * 
 *   Luke Sheneman
 *   sheneman@cs.uidaho.edu
 *
 */



#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <float.h>

#include "distclearcut.h"
#include "dmat.h"
#include "fasta.h"
#include "cmdargs.h"
#include "common.h"
#include "clearcut.h"
#include "prng.h"


/*
 * main() - 
 *
 * The entry point to the program.
 *
 */
int clearcut_main(int argc, char *argv[]) {

  DMAT *dmat;         /* The working distance matrix */
  DMAT *dmat_backup = nullptr;/* A backup distance matrix    */
  NJ_TREE *tree;      /* The phylogenetic tree       */
  NJ_ARGS *nj_args;   /* Structure for holding command-line arguments */
  long int i;

  /* some variables for tracking time */
  struct timeval tv;
  unsigned long long startUs, endUs;
  

  /* check and parse supplied command-line arguments */
  nj_args = NJ_handle_args(argc, argv);

  if(!nj_args) {
    fprintf(stderr, "Clearcut: Error processing command-line arguments.\n");
    exit(-1);
  }

  /* for verbose reporting, print the random number seed to stdout */
  if(nj_args->verbose_flag) {
    printf("PRNG SEED: %d\n", nj_args->seed);
  }

  /* Initialize Mersenne Twister PRNG */
  init_genrand(nj_args->seed);


  switch(nj_args->input_mode) {

    /* If the input type is a distance matrix */
  case NJ_INPUT_MODE_DISTANCE:

    /* parse the distance matrix */
    dmat = NJ_parse_distance_matrix(nj_args);
    if(!dmat) {
      exit(-1);
    }

    break;

    /* If the input type is a multiple sequence alignment */
  case NJ_INPUT_MODE_ALIGNED_SEQUENCES:

    /* build a distance matrix from a multiple sequence alignment */
    dmat = NJ_build_distance_matrix(nj_args);
    if(!dmat) {
      fprintf(stderr, "Clearcut: Failed to build distance matrix from alignment.\n");
      exit(-1);
    }
    
    break;

  default:

    fprintf(stderr, "Clearcut: Could not determine how to process input\n");
    exit(-1);
  }

  /*
   * Output the computed distance matrix,
   *  if the user specified one.
   */
  if(nj_args->matrixout) {
    NJ_output_matrix(nj_args, dmat);
  }
  
  /* 
   * If we are going to generate multiple trees from
   * the same distance matrix, we need to make a backup 
   * of the original distance matrix.
   */
  if(nj_args->ntrees > 1) {
    dmat_backup = NJ_dup_dmat(dmat);
  }
  
  /* process n trees */
  for(i=0;i<nj_args->ntrees;i++) {
    
    /* 
     * If the user has specified matrix shuffling, we need 
     * to randomize the distance matrix
     */
    if(nj_args->shuffle) {
      NJ_shuffle_distance_matrix(dmat);
    }

    /* RECORD THE PRECISE TIME OF THE START OF THE NEIGHBOR-JOINING */
    gettimeofday(&tv, nullptr);
    startUs = ((unsigned long long) tv.tv_sec * 1000000ULL)
      + ((unsigned long long) tv.tv_usec);

    
    /* 
     * Invoke either the Relaxed Neighbor-Joining algorithm (default)
     * or the "traditional" Neighbor-Joining algorithm 
     */
    if(nj_args->neighbor) {
      tree = NJ_neighbor_joining(nj_args, dmat);
    } else {
      tree = NJ_relaxed_nj(nj_args, dmat);
    }
  
    if(!tree) {
      fprintf(stderr, "Clearcut: Failed to construct tree.\n");
      exit(0);
    }

    /* RECORD THE PRECISE TIME OF THE END OF THE NEIGHBOR-JOINING */
    gettimeofday(&tv, nullptr);
    endUs = ((unsigned long long) tv.tv_sec * 1000000ULL)
      + ((unsigned long long) tv.tv_usec);

    /* print the time taken to perform the neighbor join */
    if(nj_args->verbose_flag) {
      if(nj_args->neighbor) {
	fprintf(stderr, "NJ tree built in %llu.%06llu secs\n",
		(endUs - startUs) / 1000000ULL,
		(endUs - startUs) % 1000000ULL);
      } else { 
	fprintf(stderr, "RNJ tree built in %llu.%06llu secs\n",
		(endUs - startUs) / 1000000ULL,
		(endUs - startUs) % 1000000ULL);
      }
    }

    /* Output the neighbor joining tree here */
    NJ_output_tree(nj_args, tree, dmat, i);
    
    NJ_free_tree(tree);  /* Free the tree */
    NJ_free_dmat(dmat);  /* Free the working distance matrix */

    /* 
     * If we need to do another iteration, lets re-initialize 
     * our working distance matrix.
     */
    if(nj_args->ntrees > 1 && i<(nj_args->ntrees-1) ) {
      dmat = NJ_dup_dmat(dmat_backup);
    }
  }
  
  /* Free the backup distance matrix */
  if(nj_args->ntrees > 1) {
    NJ_free_dmat(dmat_backup);
  }

  /* If verbosity, describe where the tree output is */
  if(nj_args->verbose_flag) {
    if(nj_args->neighbor) {
      printf("NJ tree(s) in %s\n", nj_args->outfilename);
    } else {
      printf("Relaxed NJ tree(s) in %s\n", nj_args->outfilename);
    }
  }
  
	return 0;
}





/*
 * NJ_find_hmin() - Find minimum transformed values along horizontal
 * 
 * 
 * INPUTS:
 * -------
 *      dmat -- The distance matrix
 *         a -- The index of the specific taxon in the distance matrix
 *
 * RETURNS:
 * --------
 *   <float> -- The value of the selected minimum
 *       min -- Used to transport the index of the minima out 
 *              of the function (by reference)
 * hmincount -- Return the number of minima along the horizontal
 *              (by reference)
 *
 *
 * DESCRIPTION:
 * ------------
 *
 * A fast, inline function to find the smallest transformed value 
 * along the "horizontal" portion of an entry in a distance matrix.
 *
 * Distance matrices are stored internally as continguously-allocated
 * upper-diagonal structures.  With the exception of the taxa at
 * row 0 of this upper-diagonal matrix, all taxa have both a horizontal
 * and vertical component in the distance matrix.  This function
 * scans the horizonal portion of the entry in the distance matrix
 * for the specified taxon and finds the minimum transformed value
 * along that horizontal component.
 *
 * Since multiple minima can exist along the horizontal portion
 * of the entry, I consider all minima and break ties
 * stochastically to help avoid systematic bias.
 *
 * Just searching along the horizontal portion of a row is very fast
 * since the data is stored linearly and contiguously in memory and 
 * cache locality is exploited in the distance matrix representation. 
 *
 * Look at nj.h for more information on how the distance matrix 
 * is architected.
 * 
 */
static inline
float
NJ_find_hmin(DMAT *dmat,
	     long int a,
	     long int *min,
	     long int *hmincount) {

  long int i;     /* index variable for looping                    */
  int size;       /* current size of distance matrix               */
  int mindex = 0; /* holds the current index to the chosen minimum */
  float curval;   /* used to hold current transformed values       */
  float hmin;     /* the value of the transformed minimum          */

  float *ptr, *r2, *val;  /* pointers used to reduce dereferencing in inner loop */

  /* values used for stochastic selection among multiple minima */
  float p, x;  
  long int smallcnt;

  /* initialize the min to something large */
  hmin = (float)HUGE_VAL;

  /* setup some pointers to limit dereferencing later */
  r2       = dmat->r2;
  val      = dmat->val;
  size     = dmat->size;

  /* initialize values associated with minima tie breaking */
  p        = 1.0;
  smallcnt = 0;
  
  
  ptr = &(val[NJ_MAP(a, a+1, size)]);   /* at the start of the horiz. part */
  for(i=a+1;i<size;i++) {

    curval = *(ptr++) - (r2[a] + r2[i]);  /* compute transformed distance */
    
    if(NJ_FLT_EQ(curval, hmin)) {  /* approx. equal */
      
      smallcnt++;

      p = 1.0/(float)smallcnt;
      x = genrand_real2();
      
      /* select this minimum in a way which is proportional to 
	 the number of minima found along the row so far */
      if( x < p ) {
	mindex = i;
      }

    } else if (curval < hmin) {

      smallcnt = 1;
      hmin = curval;
      mindex = i;
    }
  }
  
  /* save off the the minimum index to be returned via reference */
  *min = mindex;
  
  /* save off the number of minima */
  *hmincount = smallcnt;
  
  /* return the value of the smallest tranformed distance */
  return(hmin);
}








/*
 * NJ_find_vmin() - Find minimum transformed distance along vertical
 *
 * 
 * INPUTS:
 * -------
 *      dmat -- The distance matrix
 *         a -- The index of the specific taxon in the distance matrix
 *
 *
 * RETURNS:
 * --------
 *   <float> -- The value of the selected minimum
 *       min -- Used to transport the index of the minima out 
 *              of the function (by reference)
 * vmincount -- The number of minima along the vertical
 *              return by reference.
 *
 * DESCRIPTION:
 * ------------
 *
 * A fast, inline function to find the smallest transformed value 
 * along the "vertical" portion of an entry in a distance matrix.
 *
 * Distance matrices are stored internally as continguously-allocated
 * upper-diagonal matrices.  With the exception of the taxa at
 * row 0 of this upper-diagonal matrix, all taxa have both a horizontal
 * and vertical component in the distance matrix.  This function
 * scans the vertical portion of the entry in the distance matrix
 * for the specified taxon and finds the minimum transformed value
 * along that vertical component.
 *
 * Since multiple minima can exist along the vertical portion
 * of the entry, I consider all minima and break ties
 * stochastically to help avoid systematic bias.
 *
 * Due to cache locality reasons, searching along the vertical
 * component is going to be considerably slower than searching
 * along the horizontal.
 *
 * Look at nj.h for more information on how the distance matrix 
 * is architected.
 * 
 */
static inline
float
NJ_find_vmin(DMAT *dmat,
	     long int a,
	     long int *min,
	     long int *vmincount) {

  long int i;         /* index variable used for looping */
  long int size;      /* track the size of the matrix    */
  long int mindex = 0;/* track the index to the minimum  */
  float curval;       /* track value of current transformed distance  */
  float vmin;         /* the index to the smallest "vertical" minimum */

  /* pointers which are used to reduce pointer dereferencing in inner loop */
  float *ptr, *r2, *val;

  /* values used in stochastically breaking ties */
  float p, x;
  long int smallcnt;

  /* initialize the vertical min to something really big */
  vmin = (float)HUGE_VAL;

  /* save off some values to limit dereferencing later */
  r2       = dmat->r2;
  val      = dmat->val;
  size     = dmat->size;

  p        = 1.0;
  smallcnt = 0;

  /* start on the first row and work down */
  ptr = &(val[NJ_MAP(0, a, size)]);  
  for(i=0;i<a;i++) {

    curval = *ptr - (r2[i] + r2[a]);  /* compute transformed distance */
    
    if(NJ_FLT_EQ(curval, vmin)) {  /* approx. equal */
      
      smallcnt++;

      p = 1.0/(float)smallcnt;
      x = genrand_real2();

      /* break ties stochastically to avoid systematic bias */
      if( x < p ) {
	mindex = i;
      }

    } else if (curval < vmin) {

      smallcnt = 1;
      vmin = curval;
      mindex = i;
    }

    /* increment our working pointer to the next row down */
    ptr += size-i-1;
  }

  /* pass back the index to the minimum found so far (by reference) */
  *min = mindex;
  
  /* pass back the number of minima along the vertical */
  *vmincount = smallcnt;

  /* return the value of the smallest transformed distance */
  return(vmin);
}




/*
 * NJ_permute() - Generate random permutation using the provably unbiased
 *                Fisher-Yates Shuffle.
 *
 * INPUTS:
 * -------
 *   perm -- A pointer to the array of long ints which will be filled.
 *   size -- the length of the permutation vector 
 *
 *
 * OUTPUTS:
 * -------
 *   NONE
 *
 *
 * DESCRIPTION:
 * ------------
 *
 * Return a permuted list of numbers from 0 through size.
 * This is accomplished by initializing the permutation 
 * as an ordered list of integers and then iterating 
 * through and swapping two values selected according to the
 * Fisher-Yates method.
 *
 * This unbiased method for random permutation generation is 
 * discussed in:
 *
 *     Donald E. Knuth, The Art of Computer Programming, 
 *     Addison-Wesley, Volumes 1, 2, and 3, 3rd edition, 1998
 *
 */
static inline
void
NJ_permute(long int *perm,
	   long int size) {
  
  long int i;     /* index used for looping */
  long int swap;  /* we swap values to generate permutation */
  long int tmp;   /* used for swapping values */


  /* check to see if vector of long ints is valid */
  if(!perm) {
    fprintf(stderr, "Clearcut: nullptr permutation pointer in NJ_permute()\n");
    exit(-1);
  }
  
  /* init permutation as an ordered list of integers */
  for(i=0;i<size;i++) {
    perm[i] = i;
  }

  /* 
   * Iterate across the array from i = 0 to size -1, swapping ith element
   * with a randomly chosen element from a changing range of possible values
   */
  for(i=0;i<size;i++) {

    /* choose which element we will swap with */
    swap = i + NJ_genrand_int31_top(size-i);

    /* swap elements here */
    if(i != swap) {
      tmp        = perm[swap];
      perm[swap] = perm[i];
      perm[i]    = tmp;
    }
  }
  
  return;
}





/*
 * NJ_compute_r() - Compute post-join changes to r-vector.  In this
 *                  case, we decrement all of the accumulated distances
 *                  in the r-vector for the two nodes that were
 *                  recently joined (i.e.  x, y)
 *
 * INPUTS:
 * -------
 *   dmat -- The distance matrix 
 *      a -- The index of one of the taxa that were joined
 *      b -- The index of the other taxa that was joined
 *
 * RETURNS:
 * --------
 *   NONE
 * 
 * DESCRIPTION:
 * ------------
 *   
 * This vector of floats is used as a summary of overall distances from
 * each entry in the distance matrix to every other entry.  These values
 * are then used when computing the transformed distances from which
 * decisions concerning joining are made.
 *
 * For speed, we don't recompute r from scratch.  Instead, we decrement
 * all entries in r by the appropriate amount.  That is, r[i] -= dist(i, a)
 * and r[i] -= dist(i, b).
 *
 * As a speed optimization, I process the rows altogether for cache locality
 * purposes, and then process columns.
 *
 * The processing of the scaled r matrix (r2) is handled on-the-fly elsewhere.
 *  
 */
static inline
void
NJ_compute_r(DMAT *dmat,
	     long int a,
	     long int b) {
  
  long int i;         /* a variable used in indexing */
  float *ptrx, *ptry; /* pointers into the distance matrix */
  
  /* some variables to limit pointer dereferencing in loop */
  long int size; 
  float *r, *val;
  
  /* to limit pointer dereferencing */
  size = dmat->size;
  val  = dmat->val;
  r    = dmat->r+a+1;

  /* 
   * Loop through the rows and decrement the stored r values 
   * by the distances stored in the rows and columns of the distance 
   * matrix which are being removed post-join.
   *
   * We do the rows altogether in order to benefit from cache locality.
   */
  ptrx = &(val[NJ_MAP(a, a+1, size)]); 
  ptry = &(val[NJ_MAP(b, b+1, size)]); 

  for(i=a+1;i<size;i++) {
    *r -= *(ptrx++);  

    if(i>b) {
      *r -= *(ptry++); 
    }

    r++;
  }

  /* Similar to the above loop, we now do the columns */
  ptrx = &(val[NJ_MAP(0, a, size)]);  
  ptry = &(val[NJ_MAP(0, b, size)]);  
  r = dmat->r;
  for(i=0;i<b;i++) {
    if(i<a) {
      *r -= *ptrx;
      ptrx += size-i-1;
    }

    *r -= *ptry;
    ptry += size-i-1;
    r++;
  }

  return;
}





/*
 * NJ_check_additivity() - Check to make sure that addivity preserved by join
 *
 *
 * INPUTS:
 * -------
 *    dmat -- distance matrix
 *       a -- index into dmat for one of the rows to be joined
 *       b -- index into dmat for another row to be joined
 *
 * OUTPUTS:
 * --------
 *     int    1 if join adheres to additivity constraint
 *            0 if join does breaks additivity
 *
 * DESCRIPTION:
 * ------------
 *
 * Here we perform the check to make sure that by joining a and b we do not 
 * also break consistency (i.e. preserves additivity) with the distances between 
 * the taxa in the new clade and other nodes in the tree.  This is done quite
 * efficiently by looking up the untransformed distance between node b and 
 * some other "target" taxa in the distance matrix (which is not a nor b) and 
 * comparing that distance to the distance computed by finding the distance 
 * from node a to the proposed internal node "x" which joins (a,b).  
 *
 * If dist(x,b) + dist (b, target) == dist(b, target) then additivity is 
 * preserved, otherwise, additivity is not preserved.  If we are in 
 * additivity mode, this join should be rejected.
 *
 */
static inline
int
NJ_check_additivity(DMAT *dmat,
		    long int a,
		    long int b) {

  float a2clade, b2clade;
  float clade_dist;
  long int target;


  /* determine target taxon here */
  if(b == dmat->size-1) {
    /* if we can't do a row here, lets do a column */
    if(a==0) {
      if(b==1) {
	target = 2;
      } else {
	target = 1;
      }
    } else {
      target = 0;
    }
  } else {
    target = b+1;
  }


  /* distance between a and the root of clade (a,b) */
  a2clade = 
    ( (dmat->val[NJ_MAP(a, b, dmat->size)]) + 
      (dmat->r2[a] - dmat->r2[b]) ) / 2.0;  
  
  /* distance between b and the root of clade (a,b) */
  b2clade = 
    ( (dmat->val[NJ_MAP(a, b, dmat->size)]) + 
      (dmat->r2[b] - dmat->r2[a]) ) / 2.0;  

  /* distance between the clade (a,b) and the target taxon */
  if(b<target) {

    /* compute the distance from the clade root to the target */
    clade_dist = 
      ( (dmat->val[NJ_MAP(a, target, dmat->size)] - a2clade) +
	(dmat->val[NJ_MAP(b, target, dmat->size)] - b2clade) ) / 2.0;
    
    /* 
     * Check to see that distance from clade root to target + distance from 
     *  b to clade root are equal to the distance from b to the target 
     */
    if(NJ_FLT_EQ(dmat->val[NJ_MAP(b, target, dmat->size)], 
		 (clade_dist + b2clade))) {
      return(1);  /* join is legitimate   */
    } else {
      return(0);  /* join is illigitimate */
    }
    
  } else {

    /* compute the distance from the clade root to the target */
    clade_dist = 
      ( (dmat->val[NJ_MAP(target, a, dmat->size)] - a2clade) +
	(dmat->val[NJ_MAP(target, b, dmat->size)] - b2clade) ) / 2.0;

    /* 
     * Check to see that distance from clade root to target + distance from 
     *  b to clade root are equal to the distance from b to the target 
     */
    if(NJ_FLT_EQ(dmat->val[NJ_MAP(target, b, dmat->size)], 
		 (clade_dist + b2clade))) {
      return(1);  /* join is legitimate   */
    } else {
      return(0);  /* join is illegitimate */
    }
  }
}







/*
 * NJ_check() - Check to see if two taxa can be joined

 *
 * INPUTS:
 * -------
 * nj_args    -- Pointer to the data structure holding command-line args
 *    dmat    -- distance matrix
 *       a    -- index into dmat for one of the rows to be joined
 *       b    -- index into dmat for another row to be joined
 *     min    -- the minimum value found
 * additivity -- a flag (0 = not additive mode, 1 = additive mode)
 *
 * OUTPUTS:
 * --------
 *     int    1 if join is okay
 *            0 if join is not okay
 *
 * DESCRIPTION:
 * ------------
 *
 * This function ultimately takes two rows and makes sure that the 
 * intersection of those two rows, which has a transformed distance of
 * "min", is actually the smallest (or equal to the smallest) 
 * transformed distance for both rows (a, b).  If so, it returns
 * 1, else it returns 0.  
 *
 * Basically, we want to join two rows only if the minimum 
 * transformed distance on either row is at the intersection of
 * those two rows.
 *
 */
static inline
int 
NJ_check(NJ_ARGS *nj_args,
	 DMAT *dmat,
	 long int a,
	 long int b,
	 float min,
	 int additivity) {


  long int i, size;
  float *ptr, *val, *r2;
  

  /* some aliases for speed and readability reasons */
  val  = dmat->val;
  r2   = dmat->r2;
  size = dmat->size;


  /* now determine if joining a, b will result in broken distances */
  if(additivity) {
    if(!NJ_check_additivity(dmat, a, b)) {
      return(0);
    }
  }

  /* scan the horizontal of row b, punt if anything < min */
  ptr = &(val[NJ_MAP(b, b+1, size)]);  
  for(i=b+1;i<size;i++) {
    if( NJ_FLT_LT( (*ptr - (r2[b] + r2[i])), min) ) {
      return(0);
    }
    ptr++;
  }

  /* scan the vertical component of row a, punt if anything < min */
  if(nj_args->norandom) {  /* if we are doing random joins, we checked this */
    ptr = val + a;
    for(i=0;i<a;i++) {
      if( NJ_FLT_LT( (*ptr - (r2[i] + r2[a])), min) ) {
	return(0);
      }
      ptr += size-i-1;
    }
  }

  /* scan the vertical component of row b, punt if anything < min */
  ptr = val + b;
  for(i=0;i<b;i++) {
    if( NJ_FLT_LT( (*ptr - (r2[i] + r2[b])), min) && i!=a) {
      return(0);
    }
    ptr += size-i-1;
  }

  return(1);
}







/*
 * NJ_collapse() - Collapse the distance matrix by removing 
 *                 rows a and b from the distance matrix and
 *                 replacing them with a single new row which 
 *                 represents the internal node joining a and b
 *
 *
 * INPUTS:
 * -------
 *   dmat -- A pointer to the distance matrix
 * vertex -- A pointer to the vertex vector (vector of tree nodes)
 *           which is used in constructing the tree
 *      a -- An index to a row in the distance matrix from which we
 *           joined.  This row will be collapsed.
 *      b -- An index to a row in the distance matrix from which we
 *           joined.  This row will be collapsed.
 *
 * RETURNS:
 * --------
 *   NONE
 *
 *
 * DESCRIPTION:
 * ------------
 *
 * This function collapses the distance matrix in a way which optimizes
 * cache locality and ultimately gives us a speed improvement due to
 * cache.   At this point, we've decided to join rows a and b from
 * the distance matrix.  We will remove rows a and b from the distance  
 * matrix and replace them with a new row which represents the internal
 * node which joins rows a and b together. 
 * 
 * We always keep the matrix as compact as possible in order to 
 * get good performance from our cache in subsequent operations.  Cache
 * is the key to good performance here.  
 * 
 * Key Steps:
 * ----------
 * 
 *  1)  Fill the "a" row with the new distances of the internal node
 *      joining a and b to all other rows.  
 *  2)  Copy row 0 into what was row b
 *  3)  Increment the pointer to the start of the distance matrix
 *      by one row.
 *  4)  Decrement the size of the matrix by one row.
 *  5)  Do roughly the same thing to the r vector in order to
 *      keep it in sync with the distance matrix.
 *  6)  Compute the scaled r vector (r2) based on the updated
 *      r vector
 *
 * This keeps the distance matrix as compact as possible in memory, and
 * is a relatively fast operation. 
 *
 * This function requires that a < b
 *
 */
static inline
void
NJ_collapse(DMAT *dmat,
	    NJ_VERTEX *vertex,
	    long int a,
	    long int b) {


  long int i;     /* index used for looping */
  long int size;  /* size of dmat --> reduce pointer dereferencing */
  float a2clade;  /* distance from a to the new node that joins a and b */
  float b2clade;  /* distance from b to the new node that joins a and b */
  float cval;     /* stores distance information during loop */
  float *vptr;    /* pointer to elements in first row of dist matrix */
  float *ptra;    /* pointer to elements in row a of distance matrix */
  float *ptrb;    /* pointer to elements in row b of distance matrix */

  float *val, *r, *r2;  /* simply used to limit pointer dereferencing */


  /* We must assume that a < b */
  if(a >= b) {
    fprintf(stderr, "Clearcut: (a<b) constraint check failed in NJ_collapse()\n");
    exit(0);
  }

  /* some shortcuts to help limit dereferencing */
  val  = dmat->val;
  r    = dmat->r;
  r2   = dmat->r2;
  size = dmat->size;

  /* compute the distance from the clade components (a, b) to the new node */
  a2clade = 
    ( (val[NJ_MAP(a, b, size)]) + (dmat->r2[a] - dmat->r2[b]) ) / 2.0;  
  b2clade = 
    ( (val[NJ_MAP(a, b, size)]) + (dmat->r2[b] - dmat->r2[a]) ) / 2.0;  


  r[a] = 0.0;  /* we are removing row a, so clear dist. in r */

  /* 
   * Fill the horizontal part of the "a" row and finish computing r and r2 
   * we handle the horizontal component first to maximize cache locality
   */
  ptra = &(val[NJ_MAP(a,   a+1, size)]);   /* start ptra at the horiz. of a  */
  ptrb = &(val[NJ_MAP(a+1, b,   size)]);   /* start ptrb at comparable place */
  for(i=a+1;i<size;i++) {

    /* 
     * Compute distance from new internal node to others in 
     * the distance matrix.
     */
    cval = 
      ( (*ptra - a2clade) +
	(*ptrb - b2clade) ) / 2.0;

    /* incr.  row b pointer differently depending on where i is in loop */
    if(i<b) {
      ptrb += size-i-1;  /* traverse vertically  by incrementing by row */
    } else {
      ptrb++;            /* traverse horiz. by incrementing by column   */
    }

    /* assign the newly computed distance and increment a ptr by a column */
    *(ptra++) = cval;  

    /* accumulate the distance onto the r vector */
    r[a] += cval;
    r[i] += cval;
    
    /* scale r2 on the fly here */
    r2[i] = r[i]/(float)(size-3);
  }

  /* fill the vertical part of the "a" column and finish computing r and r2 */ 
  ptra = val + a;  /* start at the top of the columb for "a" */
  ptrb = val + b;  /* start at the top of the columb for "b" */
  for(i=0;i<a;i++) {

    /* 
     * Compute distance from new internal node to others in 
     * the distance matrix.
     */
    cval = 
      ( (*ptra - a2clade) + 
	(*ptrb - b2clade) ) / 2.0;
    
    /* assign the newly computed distance and increment a ptr by a column */
    *ptra = cval;

    /* accumulate the distance onto the r vector */
    r[a] += cval;
    r[i] += cval;

    /* scale r2 on the fly here */
    r2[i] = r[i]/(float)(size-3);

    /* here, always increment by an entire row */
    ptra += size-i-1;
    ptrb += size-i-1;
  }


  /* scale r2 on the fly here */
  r2[a] = r[a]/(float)(size-3);



  /* 
   * Copy row 0 into row b.  Again, the code is structured into two 
   * loops to maximize cache locality for writes along the horizontal
   * component of row b.
   */
  vptr = val;
  ptrb = val + b;
  for(i=0;i<b;i++) {
    *ptrb = *(vptr++);
    ptrb += size-i-1;
  }
  vptr++;  /* skip over the diagonal */
  ptrb = &(val[NJ_MAP(b, b+1, size)]); 
  for(i=b+1;i<size;i++) {
    *(ptrb++) = *(vptr++);
  }

  /* 
   * Collapse r here by copying contents of r[0] into r[b] and
   * incrementing pointer to the beginning of r by one row
   */
  r[b]    = r[0];
  dmat->r = r+1;


  /* 
   * Collapse r2 here by copying contents of r2[0] into r2[b] and
   * incrementing pointer to the beginning of r2 by one row
   */
  r2[b]    = r2[0];
  dmat->r2 = r2+1;

  /* increment dmat pointer to next row */
  dmat->val += size;
  
  /* decrement the total size of the distance matrix by one row */
  dmat->size--;

  return;
}









/*
 * NJ_neighbor_joining() - Perform a traditional Neighbor-Joining
 *
 *
 * INPUTS:
 * -------
 *  nj_args -- A pointer to a structure containing the command-line arguments
 *     dmat -- A pointer to the distance matrix
 *
 * RETURNS:
 * --------
 *   NJ_TREE * -- A pointer to the Neighbor-Joining tree.
 *
 * DESCRIPTION:
 * ------------
 *
 * This function performs a traditional Neighbor-Joining operation in which
 * the distance matrix is exhaustively searched for the global minimum 
 * transformed distance.  The two nodes which intersect at the global
 * minimum transformed distance are then joined and the distance
 * matrix is collapsed.  This process continues until there are only
 * two nodes left, at which point those nodes are joined.
 *
 */
NJ_TREE *
NJ_neighbor_joining(NJ_ARGS *nj_args,
		    DMAT *dmat) {

  
  NJ_TREE   *tree = nullptr;
  NJ_VERTEX *vertex = nullptr;

  long int a, b;
  float min;
    

  /* initialize the r and r2 vectors */
  NJ_init_r(dmat);

  /* allocate and initialize our vertex vector used for tree construction */
  vertex = NJ_init_vertex(dmat);
  if(!vertex) {
    fprintf(stderr, "Clearcut:  Could not initialize vertex in NJ_neighbor_joining()\n");
    return(nullptr);
  }
  
  /* we iterate until the working distance matrix has only 2 entries */
  while(vertex->nactive > 2) {
 
    /* 
     * Find the global minimum transformed distance from the distance matrix
     */
    min = NJ_min_transform(dmat, &a, &b);

    /* 
     * Build the tree by removing nodes a and b from the vertex array
     * and inserting a new internal node which joins a and b.  Collapse
     * the vertex array similarly to how the distance matrix and r and r2 
     * are compacted. 
     */
    NJ_decompose(dmat, vertex, a, b, 0);

    /* decrement the r and r2 vectors by the distances corresponding to a, b */
    NJ_compute_r(dmat, a, b);

    /* compact the distance matrix and the r and r2 vectors */
    NJ_collapse(dmat, vertex, a, b);
  }
  
  /* Properly join the last two nodes on the vertex list */
  tree = NJ_decompose(dmat, vertex, 0, 1, NJ_LAST);

  /* return the computed tree to the calling function */
  return(tree);
}








/*
 * NJ_relaxed_nj() -  Construct a tree using the Relaxed Neighbor-Joining
 *
 * INPUTS:
 * -------
 *   nj_args -- A pointer to a data structure containing the command-line args
 *      dmat -- A pointer to the distance matrix
 *
 * RETURNS:
 * --------
 *
 *   NJ_TREE * -- A pointer to a Relaxed Neighbor-Joining tree
 *
 * DESCRIPTION:
 * ------------
 *
 * This function implements the Relaxed Neighbor-Joining algorithm of
 *  Evans, J., Sheneman, L., and Foster, J. 
 *
 * Relaxed Neighbor-Joining works by choosing a local minimum transformed
 * distance when determining when to join two nodes.  (Traditional 
 * Neighbor-Joining chooses a global minimum transformed distance).
 *
 * The algorithm shares the property with traditional NJ that if the 
 * input distances are additive (self-consistent), then the algorithm
 * will manage to construct the true tree consistent with the additive
 * distances.  Additivity state is tracked and every proposed join is checked
 * to make sure it maintains additivity constraints.  If no 
 * additivity-preserving join is possible in a single pass, then the distance 
 * matrix is non-additive, and additivity checking is abandoned.  
 *
 * The algorithm will either attempt joins randomly, or it will perform joins
 * in a particular order.  The default behavior is to perform joins randomly,
 * but this can be switched off with a command-line switch.
 *
 * For randomized joins, all attempts are made to alleviate systematic bias
 * for the choice of rows to joins.  All tie breaking is done in a way which
 * is virtually free of bias.
 *
 * To perform randomized joins, a random permutation is constructed which 
 * specifies the order in which to attempt joins.  I iterate through the 
 * random permutation, and for each row in the random permutation, I find
 * the minimum transformed distance for that row.  If there are multiple 
 * minima, I break ties evenly.  For the row which intersects our 
 * randomly chosen row at the chosen minimum, if we are are still in 
 * additivity mode, I check to see if joining the two rows will break
 * our additivity constraints.  If not, I check to see if there exists 
 * a transformed distance which is smaller than the minimum found on the 
 * original row.  If there is, then we proceed through the random permutation
 * trying additional rows in the random order specified in the permutation.
 * If there is no smaller minimum transformed distance on either of the
 * two rows, then we join them, collapse the distance matrix, and compute
 * a new random permutation. 
 *
 * If the entire random permutation is traversed and no joins are possible
 * due to additivity constraints, then the distance matrix is not
 * additive, and additivity constraint-checking is disabled.
 *
 */
NJ_TREE *
NJ_relaxed_nj(NJ_ARGS *nj_args,
	      DMAT *dmat) {

  
  NJ_TREE *tree;
  NJ_VERTEX *vertex;
  long int a, b, t, bh, bv, i;
  float hmin, vmin, hvmin;
  float p, q, x;
  int join_flag;
  int additivity_mode;
  long int hmincount, vmincount;
  long int *permutation = nullptr;



  /* initialize the r and r2 vectors */
  NJ_init_r(dmat);

  additivity_mode = 1;

  /* allocate the permutation vector, if we are in randomize mode */
  if(!nj_args->norandom) {
    permutation = (long int *)calloc(dmat->size, sizeof(long int));
    if(!permutation) {
      fprintf(stderr, "Clearcut:  Memory allocation error in NJ_relaxed_nj()\n");
      return(nullptr);
    }
  }

  /* allocate and initialize our vertex vector used for tree construction */
  vertex = NJ_init_vertex(dmat);
  
  /* loop until there are only 2 nodes left to join */
  while(vertex->nactive > 2) {

    switch(nj_args->norandom) {

      /* RANDOMIZED JOINS */
    case 0:

      join_flag = 0;

      NJ_permute(permutation, dmat->size-1);
      for(i=0;i<dmat->size-1 && (vertex->nactive>2) ;i++) {

	a = permutation[i];

	/* find min trans dist along horiz. of row a */
	hmin = NJ_find_hmin(dmat, a, &bh, &hmincount);   
	if(a) {
	  /* find min trans dist along vert. of row a */
	  vmin = NJ_find_vmin(dmat, a, &bv, &vmincount); 
	} else {
	  vmin = hmin;
	  bv = bh;
	  vmincount = 0;
	}
	
	if(NJ_FLT_EQ(hmin, vmin)) {

	  /* 
	   * The minima along the vertical and horizontal are 
	   * the same.  Compute the proportion of minima along
	   * the horizonal (p) and the proportion of minima 
	   * along the vertical (q).
	   * 
	   * If the same minima exist along the horizonal and
	   * vertical, we break the tie in a way which is
	   * non-biased.  That is, we break the tie based on the
	   * proportion of horiz. minima versus vertical minima.
	   * 
	   */
	  p = (float)hmincount / ((float)hmincount + (float)vmincount);
	  q = 1.0 - p;
	  x = genrand_real2();
	  
	  if(x < p) {
	    hvmin = hmin;
	    b     = bh;
	  } else {
	    hvmin = vmin;
	    b     = bv;
	  }
	} else if(NJ_FLT_LT(hmin, vmin) ) {
	  hvmin   = hmin;
	  b       = bh;
	} else {
	  hvmin   = vmin;
	  b       = bv;
	}
	
	if(NJ_check(nj_args, dmat, a, b, hvmin, additivity_mode)) {

	  /* swap a and b, if necessary, to make sure a < b */
	  if(b < a) {
	    t = a;
	    a = b;
	    b = t;
	  }

	  join_flag = 1;
	
	  /* join taxa from rows a and b */
	  NJ_decompose(dmat, vertex, a, b, 0);

	  /* collapse matrix */
	  NJ_compute_r(dmat, a, b);
	  NJ_collapse(dmat, vertex, a, b);
	  
	  NJ_permute(permutation, dmat->size-1);
	}
      }
      
      /* turn off additivity if go through an entire cycle without joining */
      if(!join_flag) {
	additivity_mode = 0;
      }
      
      break;



      /* DETERMINISTIC JOINS */
    case 1:
      
      join_flag = 0;

      for(a=0;a<dmat->size-1 && (vertex->nactive > 2) ;) {
      
	/* find the min along the horizontal of row a */
	hmin = NJ_find_hmin(dmat, a, &b, &hmincount);
      
	if(NJ_check(nj_args, dmat, a, b, hmin, additivity_mode)) {
	
	  join_flag = 1;
	
	  /* join taxa from rows a and b */
	  NJ_decompose(dmat, vertex, a, b, 0);

	  /* collapse matrix */
	  NJ_compute_r(dmat, a, b);
	  NJ_collapse(dmat, vertex, a, b);

	  if(a) { 
	    a--; 
	  }
	
	} else {
	  a++;
	}
      }
    
      /* turn off additivity if go through an entire cycle without joining */
      if(!join_flag) {
	additivity_mode = 0;
      }

      break;
    }
 
  }  /* WHILE */

  /* Join the last two nodes on the vertex list */
  tree = NJ_decompose(dmat, vertex, 0, 1, NJ_LAST);
  
  if(nj_args->verbose_flag) {
    if(additivity_mode) {
      printf("Tree is additive\n");
    } else {
      printf("Tree is not additive\n");
    }
  }
  
  if(vertex) {
    NJ_free_vertex(vertex);
  }
  
  if(!nj_args->norandom && permutation) {
    free(permutation);
  }
  
  return(tree);
}






/* 
 * NJ_print_distance_matrix() - 
 *
 * Print a distance matrix
 *
 */
void
NJ_print_distance_matrix(DMAT *dmat) {

  long int i, j;

  printf("ntaxa: %ld\n", dmat->ntaxa);
  printf(" size: %ld\n", dmat->size);
  
  for(i=0;i<dmat->size;i++) {

    for(j=0;j<dmat->size;j++) {
      if(j>i) {
	printf("    %0.4f", dmat->val[NJ_MAP(i, j, dmat->size)]);  
      } else {
	printf("         -");
      }
    }


    if(dmat->r && dmat->r2) {
      printf("\t\t%0.4f", dmat->r[i]);    
      printf("\t%0.4f", dmat->r2[i]);


    
      printf("\n");

      for(j=0;j<dmat->size;j++) {
	if(j>i) {
	  printf("   %0.4f", dmat->val[NJ_MAP(i, j, dmat->size)] - (dmat->r2[i] + dmat->r2[j])); 
	} else {
	  printf("          ");
	}
      }
      
      printf("\n");
    }
  }
  
  printf("\n\n");
  
  return;
}







/*
 * NJ_output_tree() - 
 * 
 * A wrapper for the function that really prints the tree,
 * basically to get a newline in there conveniently.  :-)
 *
 * Print n trees, as specified in command-args
 *  using "count" variable from 0 to (n-1)
 *
 */
void
NJ_output_tree(NJ_ARGS *nj_args,
	       NJ_TREE *tree,
	       DMAT *dmat,
	       long int count) {

  FILE *fp;

  if(nj_args->stdout_flag) {
    fp = stdout;
  } else {

    if(count == 0) {
      fp = fopen(nj_args->outfilename, "w");  /* open for writing   */
    } else {
      fp = fopen(nj_args->outfilename, "a");  /* open for appending */
    }

    if(!fp) {
      fprintf(stderr, "Clearcut: Failed to open outfile %s\n", nj_args->outfilename);
      exit(-1);
    }
  }

  NJ_output_tree2(fp, nj_args, tree, tree, dmat);
  fprintf(fp, ";\n");
  
  if(!nj_args->stdout_flag) {
    fclose(fp);
  }

  return;
}





/*
 * NJ_output_tree2() - 
 * 
 *
 */
void
NJ_output_tree2(FILE *fp,
		NJ_ARGS *nj_args,
		NJ_TREE *tree,
		NJ_TREE *root,
		DMAT *dmat) {
  
  if(!tree) {
    return;
  }
	
  if(tree->taxa_index != NJ_INTERNAL_NODE) {

    if(nj_args->expblen) {
      fprintf(fp, "%s:%e", 
	      dmat->taxaname[tree->taxa_index],
	      tree->dist);
    } else {
      fprintf(fp, "%s:%f", 
	      dmat->taxaname[tree->taxa_index],
	      tree->dist);
    }
    
  } else {
    

    if(tree->left && tree->right) {
      fprintf(fp, "(");
    }
    if(tree->left) {
      NJ_output_tree2(fp, nj_args, tree->left, root, dmat);
    }

    if(tree->left && tree->right) {
      fprintf(fp, ",");
    }
    if(tree->right) {
      NJ_output_tree2(fp, nj_args, tree->right, root, dmat);
    }

    if(tree != root->left) { 
      if(tree->left && tree->right) {
	if(tree != root) {
	  if(nj_args->expblen) {
	    fprintf(fp, "):%e", tree->dist);
	  } else {
	    fprintf(fp, "):%f", tree->dist);
	  }
	} else {
	  fprintf(fp, ")");
	}
      }
    } else {
      fprintf(fp, ")");
    }
  }

  return;
}







/*
 * NJ_init_r()
 *
 * This function computes the r column in our matrix
 *
 */
void
NJ_init_r(DMAT *dmat) {

  long int i, j, size;
  long int index;
  float *r, *r2, *val;
  long int size1;
  float size2;
  
  r     = dmat->r;
  r2    = dmat->r2;
  val   = dmat->val;
  size  = dmat->size;
  size1 = size-1;
  size2 = (float)(size-2);

  index = 0;
  for(i=0;i<size1;i++) {
    index++;
    for(j=i+1;j<size;j++) {
      r[i] += val[index];
      r[j] += val[index];
      index++;
    }

    r2[i] = r[i]/size2;
  }
  
  return;
}











/*
 * NJ_init_vertex() - 
 *
 * Construct a vertex, which we will use to construct our tree 
 * in a true bottom-up approach.  The vertex construct is 
 * basically the center node in the initial star topology.
 *
 */
NJ_VERTEX *
NJ_init_vertex(DMAT *dmat) {
  
  long int i;
  NJ_VERTEX *vertex;
  
  /* allocate the vertex here */
  vertex = (NJ_VERTEX *)calloc(1, sizeof(NJ_VERTEX));
  
  /* allocate the nodes in the vertex */
  vertex->nodes        = (NJ_TREE **)calloc(dmat->ntaxa, sizeof(NJ_TREE *));
  vertex->nodes_handle = vertex->nodes;
  
  /* initialize our size and active variables */
  vertex->nactive = dmat->ntaxa;
  vertex->size    = dmat->ntaxa;
  
  /* initialize the nodes themselves */
  for(i=0;i<dmat->ntaxa;i++) {
    
    vertex->nodes[i] = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));

    vertex->nodes[i]->left  = nullptr;
    vertex->nodes[i]->right = nullptr;
    
    vertex->nodes[i]->taxa_index = i;
  }

  return(vertex);
}





/*
 * NJ_decompose() - 
 *
 * This function decomposes the star by creating new internal nodes
 * and joining two existing tree nodes to it
 *
 */
NJ_TREE *
NJ_decompose(DMAT *dmat,
	     NJ_VERTEX *vertex,
	     long int x,
	     long int y,
	     int last_flag) {

  NJ_TREE *new_node;
  float x2clade, y2clade;

  /* compute the distance from the clade components to the new node */
  if(last_flag) {
    x2clade = 
      (dmat->val[NJ_MAP(x, y, dmat->size)]);  
  } else {
    x2clade = 
      (dmat->val[NJ_MAP(x, y, dmat->size)])/2 +   
      ((dmat->r2[x] - dmat->r2[y])/2);
  }

  vertex->nodes[x]->dist = x2clade;

  if(last_flag) {
    y2clade = 
      (dmat->val[NJ_MAP(x, y, dmat->size)]);  
  } else {
    y2clade = 
      (dmat->val[NJ_MAP(x, y, dmat->size)])/2 +  
      ((dmat->r2[y] - dmat->r2[x])/2);
  }

  vertex->nodes[y]->dist = y2clade;

  /* allocate new node to connect two sub-clades */
  new_node = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));

  new_node->left  = vertex->nodes[x];
  new_node->right = vertex->nodes[y];
  new_node->taxa_index = NJ_INTERNAL_NODE;  /* this is not a terminal node, no taxa index */
  
  if(last_flag) {
    return(new_node);
  }

  vertex->nodes[x] = new_node;
  vertex->nodes[y] = vertex->nodes[0];
  
  vertex->nodes = &(vertex->nodes[1]);
  
  vertex->nactive--;

  return(new_node);
}



/*
 * NJ_print_vertex() - 
 *
 * For debugging, print the contents of the vertex
 *
 */
void
NJ_print_vertex(NJ_VERTEX *vertex) {

  long int i;

  printf("Number of active nodes: %ld\n", vertex->nactive);

  for(i=0;i<vertex->nactive;i++) {
    printf("%ld ", vertex->nodes[i]->taxa_index);
  }
  printf("\n");

  return;
}









/*
 * NJ_print_r() - 
 *
 */
void
NJ_print_r(DMAT *dmat) {
  
  long int i;
  
  printf("\n");
  for(i=0;i<dmat->size;i++) {
    printf("r[%ld] = %0.2f\n", i, dmat->r[i]);
  }
  printf("\n");

  return;
}





/*
 * NJ_print_taxanames() -
 *
 * Print taxa names here
 *
 */
void
NJ_print_taxanames(DMAT *dmat) {
  
  long int i;
  
  printf("Number of taxa: %ld\n", dmat->ntaxa);
  
  for(i=0;i<dmat->ntaxa;i++) {
    printf("%ld) %s\n", i, dmat->taxaname[i]);
  }
  
  printf("\n");

  return;
}




/* 
 * NJ_shuffle_distance_matrix() - 
 *
 * Randomize a distance matrix here
 *
 */
void
NJ_shuffle_distance_matrix(DMAT *dmat) {

  
  long int *perm      = nullptr;
  char **tmp_taxaname = nullptr;
  float *tmp_val      = nullptr;
  long int i, j;

  
  /* alloc the random permutation and a new matrix to hold the shuffled vals */
  perm         = (long int *)calloc(dmat->size, sizeof(long int));
  tmp_taxaname = (char **)calloc(dmat->size, sizeof(char *));
  tmp_val      = (float *)calloc(NJ_NCELLS(dmat->ntaxa), sizeof(float));
  if(!tmp_taxaname || !perm || !tmp_val) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_shuffle_distance_matrix()\n");
    exit(-1);
  }

  /* compute a permutation which will describe how to shuffle the matrix */
  NJ_permute(perm, dmat->size);

  for(i=0;i<dmat->size;i++) {
    for(j=i+1;j<dmat->size;j++) {

      if(perm[j] < perm[i]) {
	tmp_val[NJ_MAP(i, j, dmat->size)] = dmat->val[NJ_MAP(perm[j], perm[i], dmat->size)];
      } else {
	tmp_val[NJ_MAP(i, j, dmat->size)] = dmat->val[NJ_MAP(perm[i], perm[j], dmat->size)];
      }

    }
    
    tmp_taxaname[i] = dmat->taxaname[perm[i]];
  }

  /* free our random permutation */
  if(perm) {
    free(perm);
  }
  
  /* free the old value matrix */
  if(dmat->val) {
    free(dmat->val);
  }

  /* re-assign the value matrix pointers */
  dmat->val = tmp_val;
  dmat->valhandle = dmat->val;
  
  /* 
   * Free our old taxaname with its particular ordering
   * and re-assign to the new.
   */
  if(dmat->taxaname) {
    free(dmat->taxaname);
  }
  dmat->taxaname = tmp_taxaname;

  return;
}



/*
 * NJ_free_tree() - 
 *
 * Free a given NJ tree
 */
void
NJ_free_tree(NJ_TREE *node) {

  if(!node) {
    return;
  }
  
  if(node->left) {
    NJ_free_tree(node->left);
  }
  
  if(node->right) {
    NJ_free_tree(node->right);
  }
  
  free(node);

  return;
}









/*
 * NJ_print_permutation()
 *
 * Print a permutation
 *
 */
void
NJ_print_permutation(long int *perm,
		     long int size) {
  
  long int i;
  
  for(i=0;i<size-1;i++) {
    printf("%ld,", perm[i]);
  }
  printf("%ld\n", perm[size-1]);
  
  return;
}




/*
 * NJ_dup_dmat() - 
 * 
 * Duplicate a distance matrix
 *
 */
DMAT *
NJ_dup_dmat(DMAT *src) {
  
  long int i;
  DMAT *dest;
  
  /* allocate the resulting distance matrix */
  dest = (DMAT *)calloc(1, sizeof(DMAT));
  if(!dest) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_dup_dmat()\n");
    goto XIT_BAD;
  }

  dest->ntaxa = src->ntaxa;
  dest->size  = src->size;
  
  /* allocate space for array of pointers to taxanames */
  dest->taxaname = (char **)calloc(dest->ntaxa, sizeof(char *));
  if(!dest->taxaname) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_dup_dmat()\n");
    goto XIT_BAD;
  }

  /* allocate space for the taxanames themselves */
  for(i=0;i<src->ntaxa;i++) {
    dest->taxaname[i] = (char *)calloc(strlen(src->taxaname[i])+1, sizeof(char));
    if(!dest->taxaname[i]) {
      fprintf(stderr, "Clearcut: Memory allocation error in NJ_dup_dmat()\n");
      goto XIT_BAD;
    }
  }
  
  /* allocate space for the distance values */
  dest->val = (float *)calloc(NJ_NCELLS(src->ntaxa), sizeof(float));
  if(!dest->val) {
    fprintf(stderr, "Clearcut: Memory allocation error in NJ_dup_dmat()\n");
    goto XIT_BAD;
  }
  
  /* allocate space for the r and r2 vectors */
  dest->r  = (float *)calloc(src->ntaxa, sizeof(float));
  dest->r2 = (float *)calloc(src->ntaxa, sizeof(float));
  
  /* copy titles */
  for(i=0;i<src->ntaxa;i++) {
    strcpy(dest->taxaname[i], src->taxaname[i]);
  }
  
  /* copy values */
  memcpy(dest->val, src->valhandle, NJ_NCELLS(src->ntaxa)*sizeof(float));
  
  /* copy r and r2 */
  memcpy(dest->r,  src->rhandle,  src->ntaxa*sizeof(float));
  memcpy(dest->r2, src->r2handle, src->ntaxa*sizeof(float));
  
  /* track some memory addresses */
  dest->valhandle = dest->val;
  dest->rhandle   = dest->r;
  dest->r2handle  = dest->r2;
  
  return(dest);
  
 XIT_BAD:
  
  /* free what we may have allocated */
  NJ_free_dmat(dest);
  
  return(nullptr);
}




/*
 * NJ_free_dmat() - 
 */
void
NJ_free_dmat(DMAT *dmat) {
  
  long int i;
  
  if(dmat) {
    
    if(dmat->taxaname) {

      for(i=0;i<dmat->ntaxa;i++) {
	if(dmat->taxaname[i]) {
	  free(dmat->taxaname[i]);
	}
      }

      free(dmat->taxaname);
    }

    if(dmat->valhandle) {
      free(dmat->valhandle);
    }

    if(dmat->rhandle) {
      free(dmat->rhandle);
    }

    if(dmat->r2handle) {
      free(dmat->r2handle);
    }

    free(dmat);
  }
  
  return;
}





/*
 * NJ_free_vertex() - 
 *
 * Free the vertex data structure 
 *
 */
void
NJ_free_vertex(NJ_VERTEX *vertex) {
  
  if(vertex) {
    if(vertex->nodes_handle) {
      free(vertex->nodes_handle);
    }
    free(vertex);
  }

  return;
}









/*
 *
 * NJ_min_transform() - Find the smallest transformed value to identify 
 *                      which nodes to join.
 *
 * INPUTS:
 * -------
 *  dmat  -- The distance matrix
 *
 * RETURNS:
 * --------
 * <float> -- The minimimum transformed distance
 *   ret_i -- The row of the smallest transformed distance (by reference)
 *   ret_j -- The col of the smallest transformed distance (by reference)
 *
 *
 * DESCRIPTION:
 * ------------
 *
 * Used only with traditional Neighbor-Joining, this function checks the entire
 * working distance matrix and identifies the smallest transformed distance.
 * This requires traversing the entire diagonal matrix, which is itself a 
 * O(N^2) operation.
 *
 */
float
NJ_min_transform(DMAT *dmat,
		 long int *ret_i,
		 long int *ret_j) {

  long int i, j;   /* indices used for looping        */
  long int tmp_i = 0;/* to limit pointer dereferencing  */
  long int tmp_j = 0;/* to limit pointer dereferencing  */
  float smallest;  /* track the smallest trans. dist  */
  float curval;    /* the current trans. dist in loop */

  float *ptr;      /* pointer into distance matrix    */
  float *r2;       /* pointer to r2 matrix for computing transformed dists */
  
  smallest = (float)HUGE_VAL;

  /* track these here to limit pointer dereferencing in inner loop */
  ptr = dmat->val;
  r2  = dmat->r2;

  /* for every row */
  for(i=0;i<dmat->size;i++) {
    ptr++;  /* skip diagonal */
    for(j=i+1;j<dmat->size;j++) {   /* for every column */

      /* find transformed distance in matrix at i, j */
      curval = *(ptr++) - (r2[i] + r2[j]);

      /* if the transformed distanance is less than the known minimum */
      if(curval < smallest) {

	smallest = curval;
	tmp_i = i;
	tmp_j = j;
      }
    }
  }
  
  /* pass back (by reference) the coords of the min. transformed distance */
  *ret_i = tmp_i;
  *ret_j = tmp_j;
  
  return(smallest);  /* return the min transformed distance */
}








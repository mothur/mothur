


/*
 * clearcut.h
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


#ifndef _INC_CLEARCUT_H_
#define _INC_CLEARCUT_H_ 1

extern "C" {

#include "common.h"
#include "cmdargs.h"

#define NJ_VERSION "1.0.9"


#define NJ_INTERNAL_NODE -1
#define NJ_LAST 101

#define NJ_INPUT_MODE_UNKNOWN             0
#define NJ_INPUT_MODE_DISTANCE            100
#define NJ_INPUT_MODE_UNALIGNED_SEQUENCES 101
#define NJ_INPUT_MODE_ALIGNED_SEQUENCES   102

#define NJ_MODEL_NONE    100
#define NJ_MODEL_JUKES   101
#define NJ_MODEL_KIMURA  102




/*
 * DMAT - Distance Matrix
 *
 * This is arguably the most important structure in the
 * program.  This is the distance matrix, and it is used 
 * by many functions throughout the application.
 *
 * The matrix is architected as a contiguously allocated
 * upper-diagonal matrix of floats which include the 
 * diagonal.  
 *
 * Example:
 *
 *      0    1    2    3    4    5
 *   0 0.0  1.0  0.3  0.2  0.1  0.3
 *   1      0.0  0.3  0.2  0.1  0.8
 *   2           0.0  0.1  0.3  0.5 
 *   3                0.0  0.2  0.1
 *   4                     0.0  0.2
 *   5                          0.0
 *
 * The distance matrix shrinks with every join operation,
 * so I track the original and working size of the matrix 
 * inside the matrix.
 *
 * One fast optimization to shrink the distance matrix
 * involves incrementing the "val" pointer.  Thus, in 
 * addition to tracking the pointer to the distances,
 * I also track the original pointer to that I can 
 * free the memory associated with the working distance
 * matrix.
 *
 * This also applies to the r and r2 vectors which are
 * used to compute the transformed distances in the 
 * matrix.
 * 
 */

typedef struct _STRUCT_DMAT {

  long int ntaxa;   /* the original size of the distance matrix */
  long int size;    /* the current/effective size of the distance matrix */

  char **taxaname;  /* a pointer to an array of taxa name strings */

  float *val;       /* the distances */
  float *valhandle; /* to track the orig. pointer to free memory */

  float *r, *r2;    /* r and r2 vectors (used to compute transformed dists) */
  float *rhandle, *r2handle;  /* track orig. pointers to free memory */

} DMAT;



/*
 * NJ_TREE - The Tree Data Structure 
 *
 *
 * The tree is represented internally as a rooted 
 * binary tree.  Each internal node has a left and a right child.
 * 
 * Additionally, I track the distance between the current node
 * and that node's parent (i.e. the branch length).  
 * 
 * Finally, I track the index of the taxa for leaf nodes.
 *
 */
typedef struct _STRUCT_NJ_TREE {
  
  struct _STRUCT_NJ_TREE *left;  /* left child  */
  struct _STRUCT_NJ_TREE *right; /* right child */
  
  float dist;  /* branch length.  i.e. dist from node to parent */
  
  long int taxa_index; /* for terminal nodes, track the taxon index */

} NJ_TREE;



/*
 * NJ_VERTEX
 *
 * This structure is used for building trees.  It is a vector 
 * which, represents the center of the star when building the RNJ/NJ
 * tree through star-decomposition.
 *
 * It contains a vector of tree (node) pointers.  These pointers
 * get joined together by a new internal node, and the new internal
 * node is placed back into the vector of nodes (which is now smaller).
 *
 * To keep this vector in sync. with the shrinking matrix, parts of
 * the vector are shuffled around, and so a pointer to the originally
 * allocated vector is stored such that it can be freed from memory
 * later.
 *
 * The original and working sizes of the vector are also tracked.
 *
 */
typedef struct _STRUCT_NJ_VERTEX {
  
  NJ_TREE **nodes;
  NJ_TREE **nodes_handle;  /* original memory handle for freeing */

  long int nactive;  /* number of active nodes in the list */
  long int size;     /* the total size of the vertex */

} NJ_VERTEX;


/* some function prototypes */
int clearcut_main(int, char**);  

/* core function for performing Relaxed Neighbor Joining */
NJ_TREE *
NJ_relaxed_nj(NJ_ARGS *nj_args, DMAT *dmat);

/* function for performing traditional Neighbor-Joining */
NJ_TREE *
NJ_neighbor_joining(NJ_ARGS *nj_args, DMAT *dmat);

/* print the distance matrix (for debugging) */
void
NJ_print_distance_matrix(DMAT *dmat);

/* output the computed tree to stdout or to the specified file */
void
NJ_output_tree(NJ_ARGS *nj_args,
	       NJ_TREE *tree,
	       DMAT *dmat,
	       long int count);

/* the recursive function for outputting trees */
void
NJ_output_tree2(FILE *fp,
		NJ_ARGS *nj_args,
		NJ_TREE *tree,
		NJ_TREE *root,
		DMAT *dmat);

/* initialize vertex */
NJ_VERTEX *
NJ_init_vertex(DMAT *dmat);

/* used to decompose the star topology and build the tree */
NJ_TREE *
NJ_decompose(DMAT *dmat,
	     NJ_VERTEX *vertex,
	     long int x, 
	     long int y,
	     int last_flag);

/* print the vertex vector (for debugging) */
void
NJ_print_vertex(NJ_VERTEX *vertex);

/* print taxa names (for debugging) */
void
NJ_print_taxanames(DMAT *dmat);

/* initialize r-vector prior to RNJ/NJ */
void
NJ_init_r(DMAT *dmat);

/* print the r-vector (for debugging) */
void
NJ_print_r(DMAT *dmat);

/* shuffle the distance matrix, usually after reading in input */
void
NJ_shuffle_distance_matrix(DMAT *dmat);

/* free memory from the tree */
void
NJ_free_tree(NJ_TREE *node);

/* print permutations (for debugging) */
void
NJ_print_permutation(long int *perm,
		     long int size);

/* duplicate a distance matrix for multiple iterations */
DMAT *
NJ_dup_dmat(DMAT *src);

/* free the distance matrix */
void
NJ_free_dmat(DMAT *dmat);

/* free the vertex vector */
void
NJ_free_vertex(NJ_VERTEX *vertex);

/* for computing the global minimum transformed distance in traditional NJ */
float
NJ_min_transform(DMAT *dmat,
		 long int *ret_i,
		 long int *ret_j);

}

#endif /* _INC_CLEARCUT_H_ */








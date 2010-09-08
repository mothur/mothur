/*
 * cmdargs.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#ifdef USE_GNU
#include <getopt.h>
#else
#include "getopt_long.h"
#endif /* USE_GNU*/


#include "clearcut.h"
#include "cmdargs.h"


/*
 * NJ_handle_args() - 
 *
 */
NJ_ARGS *
NJ_handle_args(int argc,
	       char *argv[]) {
  
  static NJ_ARGS nj_args;
  int option_index, c;
  
  optind = 0;  //neccasary to read in arguments if code is run more than once

  struct option NJ_long_options[] = {

    /* These options don't set a flag */
    {"in",        required_argument, NULL, 'i'},
    {"out",       required_argument, NULL, 'o'},
    {"seed",      required_argument, NULL, 's'},
    {"matrixout", required_argument, NULL, 'm'},
    {"ntrees",    required_argument, NULL, 'n'},

    /* These options set a flag */
    {"verbose",      no_argument, &(nj_args.verbose_flag),  1},
    {"quiet",        no_argument, &(nj_args.quiet_flag),    1},
    {"distance",     no_argument, &(nj_args.input_mode),    NJ_INPUT_MODE_DISTANCE}, 
    {"alignment",    no_argument, &(nj_args.input_mode),    NJ_INPUT_MODE_ALIGNED_SEQUENCES},
    {"help",         no_argument, &(nj_args.help),          1},
    {"version",      no_argument, &(nj_args.version),       1},
    {"norandom",     no_argument, &(nj_args.norandom),      1},
    {"shuffle",      no_argument, &(nj_args.shuffle),       1},
    {"stdin",        no_argument, &(nj_args.stdin_flag),    1},
    {"stdout",       no_argument, &(nj_args.stdout_flag),   1},
    {"dna",          no_argument, &(nj_args.dna_flag),      1},
    {"DNA",          no_argument, &(nj_args.dna_flag),      1},
    {"protein",      no_argument, &(nj_args.protein_flag),  1},
    {"neighbor",     no_argument, &(nj_args.neighbor),      1},
    {"expblen",      no_argument, &(nj_args.expblen),       1},
    {"expdist",      no_argument, &(nj_args.expdist),       1}, 

    {"jukes",        no_argument, &(nj_args.jukes_flag),    1},
    {"kimura",       no_argument, &(nj_args.kimura_flag),   1},
    
    {0, 0, 0, 0}

  };
  
  /* initializes options to their default */
  nj_args.infilename       = NULL;
  nj_args.outfilename      = NULL;
  nj_args.matrixout        = NULL;
  nj_args.seed             = time(0);
  nj_args.verbose_flag     = 0;
  nj_args.quiet_flag       = 0;
  nj_args.input_mode       = NJ_INPUT_MODE_DISTANCE;
  nj_args.help             = 0;
  nj_args.version          = 0;
  nj_args.norandom         = 0;
  nj_args.shuffle          = 0;
  nj_args.stdin_flag       = 0;
  nj_args.stdout_flag      = 0;
  nj_args.dna_flag         = 0;
  nj_args.protein_flag     = 0;
  nj_args.correction_model = NJ_MODEL_NONE;
  nj_args.jukes_flag       = 0;
  nj_args.kimura_flag      = 0;
  nj_args.neighbor         = 0;
  nj_args.ntrees           = 1;
  nj_args.expblen          = 0;
  nj_args.expdist          = 0;

  while(1) {

    c = getopt_long(argc,
		    argv,
		    "i:o:s:m:n:vqduahVSIOrDPjkNeE",
		    NJ_long_options,
		    &option_index);
    if(c == -1) {
      break;
    }
//printf("%d\t%d\n", option_index, argc);
//for (int red = 0; red < argc; red++) { printf("%s\n", argv[red]); }
    switch(c) {

    case 0:
      if(NJ_long_options[option_index].flag) {
	break;
      }

      printf("option %s", NJ_long_options[option_index].name);
      if(optarg) {
	printf(" with arg %s", optarg);
      }
      printf("\n");
      break;

    case 'i':
      nj_args.infilename = optarg;
      break;

    case 'o':
      nj_args.outfilename = optarg;
      break;

    case 's':
      nj_args.seed = atoi(optarg);
      break;

    case 'm':
      nj_args.matrixout = optarg;
      break;

    case 'n':
      nj_args.ntrees = atoi(optarg);
      break;

    case 'v':
      nj_args.verbose_flag = 1;
      break; 

    case 'q':
      nj_args.quiet_flag = 1;
      break;

    case 'd':
      nj_args.input_mode = NJ_INPUT_MODE_DISTANCE;
      break;

    case 'a':
      nj_args.input_mode = NJ_INPUT_MODE_ALIGNED_SEQUENCES;
      break;

    case 'h':
      nj_args.help = 1;
      break;

    case 'V':
      nj_args.version = 1;
      break;

    case 'S':
      nj_args.shuffle = 1;
      break;
      
    case 'I':
      nj_args.stdin_flag = 1;
      break;

    case 'O':
      nj_args.stdin_flag = 1;
      break;

    case 'r':
      nj_args.norandom = 1;
      break;

    case 'D':
      nj_args.dna_flag = 1;
      break;

    case 'P':
      nj_args.protein_flag = 1;
      break;

    case 'j':
      nj_args.jukes_flag = 1;
      break;

    case 'k':
      nj_args.kimura_flag = 1;
      break;

    case 'N':
      nj_args.neighbor = 1;
      break;

    case 'e':
      nj_args.expblen = 1;
      break;

    case 'E':
      nj_args.expdist = 1;
      break;

   default:
     NJ_usage();
     exit(-1);
    }
  }
 
  if(optind < argc) {
    fprintf(stderr, "Clearcut: Unknown command-line argument:\n  --> %s\n", argv[optind]);
    NJ_usage();
    exit(-1);
  }
  
  if(nj_args.version) {
    printf("Clearcut Version: %s\n", NJ_VERSION);
    exit(0);
  }
  
  if(nj_args.help) {
    NJ_usage();
    exit(0);
  }
  
  /* if stdin & explicit filename are specified for input */
  if(nj_args.stdin_flag) {
    if(nj_args.infilename) {
      fprintf(stderr, "Clearcut:  Ambiguous input source specified.  Specify input filename OR stdin.\n");
      NJ_usage();
      exit(-1);
    }
  }

  /* if stdout & explicit filename are specified for output */
  if(nj_args.stdout_flag) {
    if(nj_args.outfilename) {
      fprintf(stderr, "Clearcut:  Ambiguous output specified.  Specify output filename OR stdout.\n");
      NJ_usage();
      exit(-1);
    }
  }

  /* if user did not specify stdin or filename, default to stdin */
  if(!nj_args.stdin_flag) {
    if(!nj_args.infilename) {

      fprintf(stderr, "Clearcut: No input file specified.  Using stdin.\n");
      nj_args.stdin_flag = 1;
    }
  }
  
  /* if user did not specify stdout or filename, default to stdout */
  if(!nj_args.stdout_flag) {
    if(!nj_args.outfilename) {
      
      fprintf(stderr, "Clearcut: No output file specified.  Using stdout.\n");
      nj_args.stdout_flag = 1;
    }
  }
  
  /* User must specify distance matrix or alignment */
  if(nj_args.input_mode == NJ_INPUT_MODE_UNKNOWN) {
    fprintf(stderr, "Clearcut: Must specify input type (--distance | --alignment)\n");
    NJ_usage();
    exit(-1);
  }

  /* do not allow protein or DNA options for distance matrix input */
  if(nj_args.input_mode == NJ_INPUT_MODE_DISTANCE) {
    if(nj_args.dna_flag || nj_args.protein_flag) {
      fprintf(stderr, "Clearcut:  Ambiguous arguments.  (--protein | --DNA) do not apply to distance \n");
      NJ_usage();
      exit(-1);
    }
  }
  
  /* make sure different filenames were specified for input and output */
  if(!nj_args.stdin_flag && !nj_args.stdout_flag) {

    if(!strcmp(nj_args.infilename, nj_args.outfilename)) {
      fprintf(stderr, "Clearcut: Input filename and output filename must be unique.\n");
      NJ_usage();
      exit(-1);
    }
  }

  /* make sure that user specifies DNA or Protein if dealing with alignment input */
  if(nj_args.input_mode == NJ_INPUT_MODE_ALIGNED_SEQUENCES) {
    if(!nj_args.dna_flag && !nj_args.protein_flag) {
      fprintf(stderr, "Clearcut: Must specify protein or DNA for alignment input.\n");
      NJ_usage();
      exit(-1);
    }
  }

  /* make sure that user does not specify both protein and DNA when dealing with alignment input */
  if(nj_args.input_mode == NJ_INPUT_MODE_ALIGNED_SEQUENCES) {
    if(nj_args.dna_flag && nj_args.protein_flag) {
      fprintf(stderr, "Clearcut: Specifying protein and DNA sequences are mutually exclusive options\n");
      NJ_usage();
      exit(-1);
    }
  }

  /* make sure verbose and quiet were not specified together */
  if(nj_args.verbose_flag && nj_args.quiet_flag) {
    fprintf(stderr, "Clearcut: Verbose and Quiet mode are mutually exclusive.\n");
    NJ_usage();
    exit(-1);
  }
  
  /* make sure that a correction model was specified only when providing an alignment */
  if(nj_args.input_mode == NJ_INPUT_MODE_DISTANCE) {
    if(nj_args.jukes_flag || nj_args.kimura_flag) {
      fprintf(stderr, "Clearcut:  Only specify correction model for alignment input.\n");
      NJ_usage();
      exit(-1);
    }
  } else {
    if(nj_args.jukes_flag && nj_args.kimura_flag) {
      fprintf(stderr, "Clearcut:  Only specify one correction model\n");
      NJ_usage();
      exit(-1);
    } else {
      if(nj_args.jukes_flag && !nj_args.kimura_flag) {
	nj_args.correction_model = NJ_MODEL_JUKES;
      } else if(nj_args.kimura_flag && !nj_args.jukes_flag) {
	nj_args.correction_model = NJ_MODEL_KIMURA;
      } else {
	nj_args.correction_model = NJ_MODEL_NONE;  /* DEFAULT */
      }
    }
  }
  
  /* make sure that the number of output trees is reasonable */
  if(nj_args.ntrees <= 0) {
    fprintf(stderr, "Clearcut: Number of output trees must be a positive integer.\n");
    NJ_usage();
    exit(-1);
  }
  
  /* 
   * make sure that if exponential distances are specified, 
   * we are dealing with alignment input
   */
  if(nj_args.expdist && nj_args.input_mode != NJ_INPUT_MODE_ALIGNED_SEQUENCES) {
    fprintf(stderr, "Clearcut: Exponential notation for distance matrix output requires that input be an alignment\n");
    NJ_usage();
    exit(-1);
  }
  
  return(&nj_args);
}





/*
 * NJ_print_args() - 
 *
 */
void
NJ_print_args(NJ_ARGS *nj_args) {
  
  char input_mode[32];
  
  switch (nj_args->input_mode) {
  case NJ_INPUT_MODE_DISTANCE:
    sprintf(input_mode, "Distance Matrix");
    break;
  case NJ_INPUT_MODE_UNALIGNED_SEQUENCES:
    sprintf(input_mode, "Unaligned Sequences");
    break;
  case NJ_INPUT_MODE_ALIGNED_SEQUENCES:
    sprintf(input_mode, "Aligned Sequences");
    break;
  default:
    sprintf(input_mode, "UNKNOWN");
    break;
  }

  printf("\n***  Command Line Arguments ***\n");
  
  printf("Input Mode: %s\n", input_mode);
  
  if(nj_args->stdin_flag) {
    printf("Input from STDIN\n");
  } else {
    printf("Input File: %s\n", nj_args->infilename);
  }

  if(nj_args->stdout_flag) {
    printf("Output from STDOUT\n");
  } else {
    printf("Output File: %s\n", nj_args->outfilename);
  }
  
  if(nj_args->input_mode != NJ_INPUT_MODE_DISTANCE) {
    if(nj_args->aligned_flag) {
      printf("Input Sequences Aligned: YES\n");
    } else {
      printf("Input Sequences Aligned:  NO\n");
    }
  }
  
  if(nj_args->verbose_flag) {
    printf("Verbose Mode: ON\n");
  } else {
    printf("Verbose Mode: OFF\n");
  }
  
  if(nj_args->quiet_flag) {
    printf("Quiet Mode: ON\n");
  } else {
    printf("Quiet Mode: OFF\n");
  }
  
  if(nj_args->seed) {
    printf("Random Seed: %d\n", nj_args->seed);
  }
  
  printf("\n*******\n");
  
  return;
}




/*
 * NJ_usage() -
 *
 * Print a usage message
 *
 */
void
NJ_usage(void) {
  
  printf("Usage: clearcut --in=<infilename> --out=<outfilename> [options]...\n");
  printf("GENERAL OPTIONS:\n");
  printf("  -h, --help         Display this information.\n");
  printf("  -V, --version      Print the version of this program.\n");
  printf("  -v, --verbose      More output. (Default: OFF)\n");
  printf("  -q, --quiet        Silent operation. (Default: ON)\n");
  printf("  -s, --seed=<seed>  Explicitly set the PRNG seed to a specific value.\n");
  printf("  -r, --norandom     Attempt joins deterministically.  (Default: OFF)\n");
  printf("  -S, --shuffle      Randomly shuffle the distance matrix.  (Default: OFF)\n");
  printf("  -N, --neighbor     Use traditional Neighbor-Joining algorithm. (Default: OFF)\n");

  printf("\n");
  printf("INPUT OPTIONS:\n");
  printf("  -I, --stdin        Read input from STDIN.\n");
  printf("  -d, --distance     Input file is a distance matrix. (Default: ON)\n");
  printf("  -a, --alignment    Input file is a set of aligned sequences. (Default: OFF)\n");
  printf("  -D, --DNA          Input alignment are DNA sequences.\n");
  printf("  -P, --protein      Input alignment are protein sequences.\n");

  printf("\n");
  printf("CORRECTION MODEL FOR COMPUTING DISTANCE MATRIX (Default: NO Correction):\n");
  printf("  -j, --jukes        Use Jukes-Cantor correction for computing distance matrix.\n");
  printf("  -k, --kimura       Use Kimura correction for distance matrix.\n");
  
  printf("\n");
  printf("OUTPUT OPTIONS:\n");
  printf("  -O, --stdout           Output tree to STDOUT.\n");
  printf("  -m, --matrixout=<file> Output distance matrix to specified file.\n");
  printf("  -n, --ntrees=<n>       Output n trees.  (Default: 1)\n");
  printf("  -e, --expblen          Exponential notation for branch lengths. (Default: OFF)\n");
  printf("  -E, --expdist          Exponential notation in distance output. (Default: OFF)\n");
  
  printf("\n");
  printf("EXAMPLES:\n");
  printf("  Compute tree by supplying distance matrix via stdin:\n");
  printf("  clearcut --distance < distances.txt > treefile.tre\n");
  printf("\n");
  printf("  Compute tree by supplying an alignment of DNA sequences from a file:\n");
  printf("  clearcut --alignment --DNA --in=alignment.txt --out=treefile.tre\n");
  
  return;
}




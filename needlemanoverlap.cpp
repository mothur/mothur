/*
 *  needleman.cpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This class is an Alignment child class that implements the Gotoh pairwise alignment algorithm as described in:
 *		
 *		Gotoh O. 1982.  An improved algorithm for matching biological sequences.  J. Mol. Biol.  162:705-8.
 *		Myers, EW & Miller, W.  1988.  Optimal alignments in linear space.  Comput Appl Biosci. 4:11-7.
 *
 *	This method is nice because it allows for an affine gap penalty to be assessed, which is analogous to what is used
 *	in blast and is an alternative to Needleman-Wunsch, which only charges the same penalty for each gap position.
 *	Because this method typically has problems at the ends when two sequences do not full overlap, we employ a separate
 *	method to fix the ends (see Overlap class documentation)
 *
 */

#include "alignmentcell.hpp"
#include "alignment.hpp"
#include "overlap.hpp"
#include "needlemanoverlap.hpp"

/**************************************************************************************************/

NeedlemanOverlap::NeedlemanOverlap(float gO, float m, float mm, int r) ://	note that we don't have a gap extend
gap(gO), match(m), mismatch(mm), Alignment(r) {							//	the gap openning penalty is assessed for
																		//	every gapped position
	for(int i=1;i<nCols;i++){
		alignment[0][i].prevCell = 'l';					//	initialize first row by pointing all poiters to the left
		alignment[0][i].cValue = 0;						//	and the score to zero
	}
	
	for(int i=1;i<nRows;i++){
		alignment[i][0].prevCell = 'u';					//	initialize first column by pointing all poiters upwards
		alignment[i][0].cValue = 0;						//	and the score to zero
	}
}

/**************************************************************************************************/

NeedlemanOverlap::~NeedlemanOverlap(){	/*	do nothing	*/	}

/**************************************************************************************************/

void NeedlemanOverlap::align(string A, string B){
	
	seqA = ' ' + A;	lA = seqA.length();		//	algorithm requires a dummy space at the beginning of each string
	seqB = ' ' + B;	lB = seqB.length();		//	algorithm requires a dummy space at the beginning of each string
	
	for(int i=1;i<lB;i++){					//	This code was largely translated from Perl code provided in Ex 3.1 
		for(int j=1;j<lA;j++){				//	of the O'Reilly BLAST book.  I found that the example output had a
											//	number of errors
			float diagonal;
			if(seqB[i] == seqA[j])	{	diagonal = alignment[i-1][j-1].cValue + match;		}
			else					{	diagonal = alignment[i-1][j-1].cValue + mismatch;	}
			
			float up	= alignment[i-1][j].cValue + gap;
			float left	= alignment[i][j-1].cValue + gap;
			
			if(diagonal >= up){
				if(diagonal >= left){
					alignment[i][j].cValue = diagonal;
					alignment[i][j].prevCell = 'd';
				}
				else{
					alignment[i][j].cValue = left;
					alignment[i][j].prevCell = 'l';
				}
			}
			else{
				if(up >= left){
					alignment[i][j].cValue = up;
					alignment[i][j].prevCell = 'u';
				}
				else{
					alignment[i][j].cValue = left;
					alignment[i][j].prevCell = 'l';
				}
			}
		}
	}
	Overlap over;						
	over.setOverlap(alignment, lA, lB, 0);		//	Fix gaps at the beginning and end of the sequences
	traceBack();								//	Traceback the alignment to populate seqAaln and seqBaln
}

/**************************************************************************************************/


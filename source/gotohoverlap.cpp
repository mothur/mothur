/*
 *  gotohoverlap.cpp
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
#include "overlap.hpp"
#include "alignment.hpp"
#include "gotohoverlap.hpp"

/**************************************************************************************************/

GotohOverlap::GotohOverlap(float gO, float gE, float f, float mm, int r) :
	gapOpen(gO), gapExtend(gE), match(f), mismatch(mm), Alignment(r) {
	
	try {
		for(int i=1;i<nCols;i++){				//	we initialize the dynamic programming matrix by setting the pointers in
			alignment[0][i].prevCell = 'l';		//	the first row to the left
			alignment[0][i].cValue = 0;
			alignment[0][i].dValue = 0;
		}
		
		for(int i=1;i<nRows;i++){				//	we initialize the dynamic programming matrix by setting the pointers in
			alignment[i][0].prevCell = 'u';		//	the first column upward
			alignment[i][0].cValue = 0;
			alignment[i][0].iValue = 0;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "GotohOverlap", "GotohOverlap");
		exit(1);
	}
}

/**************************************************************************************************/

void GotohOverlap::align(string A, string B){
	try {
		seqA = ' ' + A;	lA = seqA.length();		//	the algorithm requires that the first character be a dummy value
		seqB = ' ' + B;	lB = seqB.length();		//	the algorithm requires that the first character be a dummy value
		
		for(int i=1;i<lB;i++){					//	the recursion here is shown in Webb and Miller, Fig. 1A.  Note that 
			for(int j=1;j<lA;j++){				//	if we need to conserve on space we should see Fig. 1B, which is linear
				//	in space, which I think is unnecessary
				float diagonal;
				if(seqB[i] == seqA[j])	{	diagonal = alignment[i-1][j-1].cValue + match;		}
				else					{	diagonal = alignment[i-1][j-1].cValue + mismatch;	}
				
				alignment[i][j].iValue = max(alignment[i][j-1].iValue, alignment[i][j-1].cValue + gapOpen) + gapExtend;
				alignment[i][j].dValue = max(alignment[i-1][j].dValue, alignment[i-1][j].cValue + gapOpen) + gapExtend;
				
				if(alignment[i][j].iValue > alignment[i][j].dValue){
					if(alignment[i][j].iValue > diagonal){
						alignment[i][j].cValue = alignment[i][j].iValue;
						alignment[i][j].prevCell = 'l';
					}
					else{
						alignment[i][j].cValue = diagonal;
						alignment[i][j].prevCell = 'd';
					}
				}
				else{
					if(alignment[i][j].dValue > diagonal){
						alignment[i][j].cValue = alignment[i][j].dValue;
						alignment[i][j].prevCell = 'u';
					}
					else{
						alignment[i][j].cValue = diagonal;
						alignment[i][j].prevCell = 'd';
					}
				}
				
			}
		}
		Overlap over;
		over.setOverlap(alignment, lA, lB, 0);	//	Fix the gaps at the ends of the sequences
		traceBack();							//	Construct the alignment and set seqAaln and seqBaln
		
	}
	catch(exception& e) {
		m->errorOut(e, "GotohOverlap", "align");
		exit(1);
	}
}

/**************************************************************************************************/

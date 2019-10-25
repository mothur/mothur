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

NeedlemanOverlap::NeedlemanOverlap(float gO, float f, float mm, int r) ://	note that we don't have a gap extend
gap(gO), match(f), mismatch(mm), Alignment(r) {							//	the gap openning penalty is assessed for
	try {																	//	every gapped position
		for(int i=1;i<nCols;i++){
			alignment[0][i].prevCell = 'l';					//	initialize first row by pointing all poiters to the left
			alignment[0][i].cValue = 0;						//	and the score to zero
		}
		
		for(int i=1;i<nRows;i++){
			alignment[i][0].prevCell = 'u';					//	initialize first column by pointing all poiters upwards
			alignment[i][0].cValue = 0;						//	and the score to zero
		}
	
	}
	catch(exception& e) {
		m->errorOut(e, "NeedlemanOverlap", "NeedlemanOverlap");
		exit(1);
	}
}
/**************************************************************************************************/

NeedlemanOverlap::~NeedlemanOverlap(){	/*	do nothing	*/	}

/**************************************************************************************************/

void NeedlemanOverlap::align(string A, string B, bool constructMaps) {
	try {
	
		seqA = ' ' + A;	lA = seqA.length();		//	algorithm requires a dummy space at the beginning of each string
		seqB = ' ' + B;	lB = seqB.length();		//	algorithm requires a dummy space at the beginning of each string

		const float gap = this->gap;
		const float match = this->match;
		const float mismatch = this->mismatch;
		const char* const seqARef = seqA.c_str();
		const char* const seqBRef = seqB.c_str();

		if (lA > nRows) { m->mothurOut("One of your candidate sequences is longer than you longest template sequence. Your longest template sequence is " + toString(nRows) + ". Your candidate is " + toString(lA) + ".\n");   }
		
		for(int i=1;i<lB;i++){					//	This code was largely translated from Perl code provided in Ex 3.1
            AlignmentCell *const alignCellsI = &alignment[i][0];
            AlignmentCell *const alignCellsUp = &alignment[i - 1][0];

			for(int j=1;j<lA;j++){				//	of the O'Reilly BLAST book.  I found that the example output had a
	
				//	number of errors
				float diagonal;
				if(seqARef[i] == seqBRef[j])	{	diagonal = alignCellsUp[j - 1].cValue + match;		}
				else					{	diagonal = alignCellsUp[j - 1].cValue + mismatch;	}
			
				const float up	= alignCellsUp[j].cValue + gap;
				const float left	= alignCellsI[j - 1].cValue + gap;
				
				if(diagonal >= up){
					if(diagonal >= left){
						alignCellsI[j].cValue = diagonal;
						alignCellsI[j].prevCell = 'd';
					}
					else{
						alignCellsI[j].cValue = left;
						alignCellsI[j].prevCell = 'l';
					}
				}
				else{
					if(up >= left){
						alignCellsI[j].cValue = up;
						alignCellsI[j].prevCell = 'u';
					}
					else{
						alignCellsI[j].cValue = left;
						alignCellsI[j].prevCell = 'l';
					}
				}
			}
		}

		Overlap over;						
		over.setOverlap(alignment, lA, lB, 0);		//	Fix gaps at the beginning and end of the sequences
		traceBack(constructMaps);					//	Traceback the alignment to populate seqAaln and seqBaln
	
	}
	catch(exception& e) {
		m->errorOut(e, "NeedlemanOverlap", "align");
		exit(1);
	}

}
/**************************************************************************************************/

void NeedlemanOverlap::alignPrimer(string A, string B, bool constructMaps) {
	try {
        
		seqA = ' ' + A;	lA = seqA.length();		//	algorithm requires a dummy space at the beginning of each string
		seqB = ' ' + B;	lB = seqB.length();		//	algorithm requires a dummy space at the beginning of each string
        
		const float gap = this->gap;
		const float match = this->match;
		const float mismatch = this->mismatch;
		const char* const seqARef = seqA.c_str();
		const char* const seqBRef = seqB.c_str();
		
		if (lA > nRows) { m->mothurOut("One of your candidate sequences is longer than you longest template sequence. Your longest template sequence is " + toString(nRows) + ". Your candidate is " + toString(lA) + "."); m->mothurOutEndLine();  }
		
		for(int i=1;i<lB;i++){					//	This code was largely translated from Perl code provided in Ex 3.1
            
			AlignmentCell *const alignCellsI = &alignment[i][0];
			AlignmentCell *const alignCellsISub1 = &alignment[i - 1][0];
			
			for(int j=1;j<lA;j++){				//	of the O'Reilly BLAST book.  I found that the example output had a
                
				//	number of errors
				float diagonal;
				if(isEquivalent(seqB[i],seqA[j])) {	diagonal = alignCellsISub1[j - 1].cValue + match; }
				else {	diagonal = alignCellsISub1[j - 1].cValue + mismatch;	}
                
				const float up	= alignCellsISub1[j].cValue + gap;
				const float left	= alignCellsI[j - 1].cValue + gap;
				
				if(diagonal >= up){
					if(diagonal >= left){
						alignCellsI[j].cValue = diagonal;
						alignCellsI[j].prevCell = 'd';
					}
					else{
						alignCellsI[j].cValue = left;
						alignCellsI[j].prevCell = 'l';
					}
				}
				else{
					if(up >= left){
						alignCellsI[j].cValue = up;
						alignCellsI[j].prevCell = 'u';
					}
					else{
						alignCellsI[j].cValue = left;
						alignCellsI[j].prevCell = 'l';
					}
				}
			}
		}
        
		Overlap over;
		over.setOverlap(alignment, lA, lB, 0);		//	Fix gaps at the beginning and end of the sequences
		traceBack(constructMaps);					//	Traceback the alignment to populate seqAaln and seqBaln
        
	}
	catch(exception& e) {
		m->errorOut(e, "NeedlemanOverlap", "alignPrimer");
		exit(1);
	}
    
}
//********************************************************************/
bool NeedlemanOverlap::isEquivalent(char oligo, char seq){
	try {
		
        bool same = true;
					
        oligo = toupper(oligo);
        seq = toupper(seq);
        
        if(oligo != seq){
            if(oligo == 'A' && (seq != 'A' && seq != 'M' && seq != 'R' && seq != 'W' && seq != 'D' && seq != 'H' && seq != 'V'))       {	same = false;	}
            else if(oligo == 'C' && (seq != 'C' && seq != 'Y' && seq != 'M' && seq != 'S' && seq != 'B' && seq != 'H' && seq != 'V'))       {	same = false;	}
            else if(oligo == 'G' && (seq != 'G' && seq != 'R' && seq != 'K' && seq != 'S' && seq != 'B' && seq != 'D' && seq != 'V'))       {	same = false;	}
            else if(oligo == 'T' && (seq != 'T' && seq != 'Y' && seq != 'K' && seq != 'W' && seq != 'B' && seq != 'D' && seq != 'H'))       {	same = false;	}
            else if((oligo == '.' || oligo == '-'))           {	same = false;	}
            else if((oligo == 'N' || oligo == 'I') && (seq == 'N'))                         {	same = false;	}
            else if(oligo == 'R' && (seq != 'A' && seq != 'G'))                        {	same = false;	}
            else if(oligo == 'Y' && (seq != 'C' && seq != 'T'))                        {	same = false;	}
            else if(oligo == 'M' && (seq != 'C' && seq != 'A'))                        {	same = false;	}
            else if(oligo == 'K' && (seq != 'T' && seq != 'G'))                        {	same = false;	}
            else if(oligo == 'W' && (seq != 'T' && seq != 'A'))                        {	same = false;	}
            else if(oligo == 'S' && (seq != 'C' && seq != 'G'))                        {	same = false;	}
            else if(oligo == 'B' && (seq != 'C' && seq != 'T' && seq != 'G'))       {	same = false;	}
            else if(oligo == 'D' && (seq != 'A' && seq != 'T' && seq != 'G'))       {	same = false;	}
            else if(oligo == 'H' && (seq != 'A' && seq != 'T' && seq != 'C'))       {	same = false;	}
            else if(oligo == 'V' && (seq != 'A' && seq != 'C' && seq != 'G'))       {	same = false;	}
        }

		
		
		return same;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimOligos", "countDiffs");
		exit(1);
	}
}
/**************************************************************************************************/


/*
 *  alignment.cpp
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *  This is a class for an abstract datatype for classes that implement various types of alignment	algorithms.
 *	As of 12/18/08 these included alignments based on blastn, needleman-wunsch, and the	Gotoh algorithms
 *  
 */

#include "alignmentcell.hpp"
#include "alignment.hpp"


/**************************************************************************************************/

Alignment::Alignment() {	/*	do nothing	*/	}

/**************************************************************************************************/

Alignment::Alignment(int A) : nCols(A), nRows(A) {
	
	alignment.resize(nRows);			//	For the Gotoh and Needleman-Wunsch we initialize the dynamic programming
	for(int i=0;i<nRows;i++){			//	matrix by initializing a matrix that is A x A.  By default we will set A
		alignment[i].resize(nCols);		//	at 2000 for 16S rRNA gene sequences
	}	
	
}

/**************************************************************************************************/

void Alignment::traceBack(){			//	This traceback routine is used by the dynamic programming algorithms
										//	to fill the values of seqAaln and seqBaln
	seqAaln = "";
	seqBaln = "";
	int row = lB-1;
	int column = lA-1;
//	seqAstart = 1;
//	seqAend = column;
	
	AlignmentCell currentCell = alignment[row][column];	//	Start the traceback from the bottom-right corner of the
														//	matrix
	
	if(currentCell.prevCell == 'x'){	seqAaln = seqBaln = "NOALIGNMENT";		}//If there's an 'x' in the bottom-
	else{												//	right corner bail out because it means nothing got aligned
		while(currentCell.prevCell != 'x'){				//	while the previous cell isn't an 'x', keep going...
			
			if(currentCell.prevCell == 'u'){			//	if the pointer to the previous cell is 'u', go up in the
				seqAaln = '-' + seqAaln;				//	matrix.  this indicates that we need to insert a gap in
				seqBaln = seqB[row] + seqBaln;			//	seqA and a base in seqB
				currentCell = alignment[--row][column];
			}
			else if(currentCell.prevCell == 'l'){		//	if the pointer to the previous cell is 'l', go to the left
				seqBaln = '-' + seqBaln;				//	in the matrix.  this indicates that we need to insert a gap
				seqAaln = seqA[column] + seqAaln;		//	in seqB and a base in seqA
				currentCell = alignment[row][--column];
			}
			else{
				seqAaln = seqA[column] + seqAaln;		//	otherwise we need to go diagonally up and to the left,
				seqBaln = seqB[row] + seqBaln;			//	here we add a base to both alignments
				currentCell = alignment[--row][--column];
			}
		}
	}
	
	pairwiseLength = seqAaln.length();
	seqAstart = 1;	seqAend = 0;
	seqBstart = 1;	seqBend = 0;
	
	for(int i=0;i<seqAaln.length();i++){
		if(seqAaln[i] != '-' && seqBaln[i] == '-')		{	seqAstart++;	}
		else if(seqAaln[i] == '-' && seqBaln[i] != '-')	{	seqBstart++;	}
		else											{	break;			}
	}
	
	pairwiseLength -= (seqAstart + seqBstart - 2);
	
	for(int i=seqAaln.length()-1; i>=0;i--){
		if(seqAaln[i] != '-' && seqBaln[i] == '-')		{	seqAend++;		}
		else if(seqAaln[i] == '-' && seqBaln[i] != '-')	{	seqBend++;		}
		else											{	break;			}
	}
	pairwiseLength -= (seqAend + seqBend);

	seqAend = seqA.length() - seqAend - 1;
	seqBend = seqB.length() - seqBend - 1;

}

/**************************************************************************************************/

string Alignment::getSeqAAln(){
	return seqAaln;										//	this is called to get the alignment of seqA
}

/**************************************************************************************************/

string Alignment::getSeqBAln(){
	return seqBaln;										//	this is called to get the alignment of seqB							
}

/**************************************************************************************************/

int Alignment::getCandidateStartPos(){
	return seqAstart;									//	this is called to report the quality of the alignment
}

/**************************************************************************************************/

int Alignment::getCandidateEndPos(){
	return seqAend;										//	this is called to report the quality of the alignment
}

/**************************************************************************************************/

int Alignment::getTemplateStartPos(){
	return seqBstart;									//	this is called to report the quality of the alignment
}

/**************************************************************************************************/

int Alignment::getTemplateEndPos(){
	return seqBend;										//	this is called to report the quality of the alignment
}

/**************************************************************************************************/

int Alignment::getPairwiseLength(){
	return pairwiseLength;								//	this is the pairwise alignment length
}

/**************************************************************************************************/

//int Alignment::getLongestTemplateGap(){
//
//	int length = seqBaln.length();
//	int longGap = 0;
//	int gapLength = 0;
//	
//	int start = seqAstart;
//	if(seqAstart < seqBstart){	start = seqBstart;	}
//	for(int i=seqAstart;i<length;i++){
//		if(seqBaln[i] == '-'){
//			gapLength++;
//		}
//		else{
//			if(gapLength > 0){
//				if(gapLength > longGap){	longGap = gapLength;	}
//			}
//			gapLength = 0;
//		}
//	}
//	return longGap;
//}

/**************************************************************************************************/

/*
 *  overlap.cpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This class cleans up the alignment at the 3' end of the alignments.  Because the Gotoh and Needleman-Wunsch
 *	algorithms start the traceback from the lower-right corner of the dynamic programming matrix, there may be a lot of
 *	scattered bases in the alignment near the 3' end of the alignment.  Here we basically look for the largest score
 *	in the last column and row to determine whether there should be exta gaps in sequence A or sequence B.  The gap
 *	issues at the 5' end of the alignment seem to take care of themselves in the traceback.
 *
 */

#include "alignmentcell.hpp"
#include "overlap.hpp"


/**************************************************************************************************/

int Overlap::maxRow(vector<vector<AlignmentCell> >& alignment, const int band){
	
	float max = -100;
	int end = lA - 1;
	int index = end;
	
	for(int i=band;i<lB;i++){					//	find the row where the right most column has the highest alignment
		if(alignment[i][end].cValue >= max){	//	score.
			index = i;
			max = alignment[i][end].cValue;
		}
	}
	return index;
}

/**************************************************************************************************/

int Overlap::maxColumn(vector<vector<AlignmentCell> >& alignment, const int band){
	
	float max = -100;
	int end = lB - 1;
	int index = end;
	
	for(int i=band;i<lA;i++){					//	find the column where the bottom most column has the highest
		if(alignment[end][i].cValue >= max){	//	alignment score.
			index = i;
			max = alignment[end][i].cValue;
		}
	}
	return index;
}

/**************************************************************************************************/

void Overlap::setOverlap(vector<vector<AlignmentCell> >& alignment, const int nA, const int nB, const int band=0){
	
	lA = nA;
	lB = nB;	
	
	int rowIndex = maxRow(alignment, band);		//	get the index for the row with the highest right hand side score
	int colIndex = maxColumn(alignment, band);	//	get the index for the column with the highest bottom row score
		
	int row = lB-1;
	int column = lA-1;
	
	if(colIndex == column && rowIndex == row){}	//	if the max values are the lower right corner, then we're good
	else if(alignment[row][colIndex].cValue < alignment[rowIndex][column].cValue){
		for(int i=rowIndex+1;i<lB;i++){			//	decide whether sequence A or B needs the gaps at the end either set 
			alignment[i][column].prevCell = 'u';//	the pointer upwards or...
		}
		
	}
	else {
		for(int i=colIndex+1;i<lA;i++){
			alignment[row][i].prevCell = 'l';	//	...to the left
		}
	}
}												//	the traceback should take care of the gaps at the 5' end

/**************************************************************************************************/


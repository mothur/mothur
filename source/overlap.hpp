#ifndef OVERLAP_H
#define OVERLAP_H

/*
 *  overlap.hpp
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

#include "mothur.h"

/**************************************************************************************************/

class Overlap {
	
public:
	Overlap(){};
	~Overlap(){};
	void setOverlap(vector<vector<AlignmentCell> >&, const int, const int, const int);
private:
	int maxRow(vector<vector<AlignmentCell> >&, const int);
	int maxColumn(vector<vector<AlignmentCell> >&, const int);
	int lA, lB;
};

/**************************************************************************************************/

#endif

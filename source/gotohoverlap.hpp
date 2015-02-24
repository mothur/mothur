#ifndef GOTOHOVERLAP_H
#define GOTOHOVERLAP_H

/*
 *  gotohoverlap.h
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

#include "mothur.h"
#include "alignment.hpp"

/**************************************************************************************************/

class GotohOverlap : public Alignment {
	
public:
	GotohOverlap(float, float, float, float, int);
	void align(string, string);
	
	~GotohOverlap() {}
	
private:
	float gapOpen;
	float gapExtend;
	float match;
	float mismatch;
};

/**************************************************************************************************/

#endif

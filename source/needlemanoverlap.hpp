#ifndef NEEDLEMAN_H
#define NEEDLEMAN_H

/*
 *  needleman.h
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This class is an Alignment child class that implements the Needleman-Wunsch pairwise alignment algorithm as
 *	described in:
 *		
 *	Needleman SB & Wunsch CD.  1970.  A general method applicable to the search for similarities in the amino acid
 *		sequence of two proteins.  J Mol Biol.  48:443-53.
 *	Korf I, Yandell M, & Bedell J.  2003.  BLAST.  O'Reilly & Associates.  Sebastopol, CA.
 *
 *	This method is simple as it assesses a consistent penalty for each gap position.  Because this method typically has
 *	problems at the ends when two sequences do not full overlap, we employ a separate method to fix the ends (see
 *	Overlap class documentation)
 *
 */

#include "mothur.h"
#include "alignment.hpp"

/**************************************************************************************************/

class NeedlemanOverlap : public Alignment {
	
public:
	NeedlemanOverlap(float, float, float, int);
	~NeedlemanOverlap();
	void align(string, string, bool createBaseMap=false);
    void alignPrimer(string, string);
	
private:	
	float gap;
	float match;
	float mismatch;
    bool isEquivalent(char, char);
};

/**************************************************************************************************/

#endif

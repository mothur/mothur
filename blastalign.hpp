/*
 *  blastalign.hpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This is a basic alignment method that gets the blast program to do the heavy lifting.  In the future, we should
 *	probably incorporate NCBI's library so that we don't have to call on a user-supplied executable.  This is a child
 *	of the Alignment class, which requires a constructor and align method.
 *
 */
 
#include "mothur.h"

class BlastAlignment : public Alignment {

public:
	BlastAlignment(float, float, float, float);
	~BlastAlignment();
	void align(string, string);
private:

	string candidateFileName;
	string templateFileName;
	string blastFileName;

	void setPairwiseSeqs();
	float match;
	float mismatch;
	float gapOpen;
	float gapExtend;
};


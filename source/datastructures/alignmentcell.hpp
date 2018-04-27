#ifndef ALIGNMENTCELL_H
#define ALIGNMENTCELL_H

/*
 *  alignmentcell.hpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This class is pretty basic.  Each AlignmentCell object contains a pointer to the previous cell and different values
 *	used to calcualte the alignment.  Initially everything is set to zero and all pointers are set to 'x'
 *  
 */
#include "mothurout.h"
//********************************************************************************************************************

class AlignmentCell {
	
public:
	AlignmentCell();
    ~AlignmentCell() {}
	char prevCell;
	float cValue;
	float dValue;
	float iValue;
};

//********************************************************************************************************************

#endif

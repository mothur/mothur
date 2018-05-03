/*
 *  alignmentcell.cpp
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This class is pretty basic.  Each AlignmentCell object contains a pointer to the previous cell and different values
 *	used to calcualte the alignment.  Initially everything is set to zero and all pointers are set to 'x'
 *  
 */

#include "alignmentcell.hpp"

//********************************************************************************************************************

AlignmentCell::AlignmentCell() : prevCell('x'), cValue(0.0000), dValue(0.0000), iValue(0.0000) {}

//********************************************************************************************************************

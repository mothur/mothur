#ifndef FORMATCOLUMN_H
#define FORMATCOLUMN_H
/*
 *  formatcolumn.h
 *  Mothur
 *
 *  Created by westcott on 1/13/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "formatmatrix.h"


/******************************************************/

class FormatColumnMatrix : public FormatMatrix {
	
public:
	FormatColumnMatrix(string);
	~FormatColumnMatrix();
	int read(NameAssignment*);
    int read(CountTable*);
	
private:
	ifstream fileHandle;
	string filename;
    Utils util;
	
};

/******************************************************/

#endif


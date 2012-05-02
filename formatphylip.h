#ifndef FORMATPHYLIP_H
#define FORMATPHYLIP_H

/*
 *  formatphylip.h
 *  Mothur
 *
 *  Created by westcott on 1/13/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "formatmatrix.h"

/******************************************************/

class FormatPhylipMatrix : public FormatMatrix {
	
public:
	FormatPhylipMatrix(string);
	~FormatPhylipMatrix();
	int read(NameAssignment*);
private:
	ifstream fileHandle;
	string filename;
};

/******************************************************/

#endif


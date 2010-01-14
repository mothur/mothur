#ifndef FORMATMATRIX_H
#define FORMATMATRIX_H

/*
 *  formatmatrix.h
 *  Mothur
 *
 *  Created by westcott on 1/13/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "listvector.hpp"
#include "nameassignment.hpp"


//**********************************************************************************************************************

class FormatMatrix {

public:
	FormatMatrix(){	}
	virtual ~FormatMatrix() {}
	
	virtual void read(NameAssignment*){};
	
	void setCutoff(float c)			{	cutoff = c;			}
	ListVector* getListVector()		{	return list;		}
	string getFormattedFileName()	{	return distFile;	}
	vector<int> getRowPositions()	{	return rowPos;		}
	
protected:
	ListVector* list;
	float cutoff;
	string distFile;
	vector<int> rowPos;
};

//**********************************************************************************************************************

#endif


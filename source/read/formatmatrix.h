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
#include "counttable.h"


//**********************************************************************************************************************
//  This class takes a distance matrix file and converts it to a file where each row contains all distances below the cutoff
//  for a given sequence.

//  Example:
	/*	5
		A  
		B  0.01 
		C  0.015 0.03 
		D  0.03 0.02 0.02  
		E  0.04 0.05 0.03 0.02   
		
		becomes 
		
		0	4	1	0.01	2	0.015	3	0.03	4	0.04	
		1	4	0	0.01	2	0.03	3	0.02	4	0.05	
		2	4	0	0.015	1	0.03	3	0.02	4	0.03	
		3	4	0	0.03	1	0.02	2	0.02	4	0.02	
		4	4	0	0.04	1	0.05	2	0.03	3	0.02
		
		column 1 - sequence name converted to row number
		column 2 - numDists under cutoff
		rest of line - sequence row -> distance, sequence row -> distance
		
		if you had a cutoff of 0.03 then the file would look like,
		
		0	3	1	0.01	2	0.015	3	0.03	
		1	3	0	0.01	2	0.03	3	0.02	
		2	4	0	0.015	1	0.03	3	0.02	4	0.03	
		3	4	0	0.03	1	0.02	2	0.02	4	0.02	
		4	2	2	0.03	3	0.02	
		
		This class also creates a vector of ints, rowPos.
		
		rowPos[0] = position in the file of distances related to sequence 0.
		If a sequence is excluded by the cutoff, it's rowPos = -1.
*/
//**********************************************************************************************************************

class FormatMatrix {

public:
	FormatMatrix(){	m = MothurOut::getInstance(); }
	virtual ~FormatMatrix() {}
	
	virtual int read(NameAssignment*){ return 1; }
    virtual int read(CountTable*){ return 1; }
	
	void setCutoff(float c)			{	cutoff = c;			}
	ListVector* getListVector()		{	return list;		}
	string getFormattedFileName()	{	return distFile;	}
	vector<int> getRowPositions()	{	return rowPos;		}
	
protected:
	ListVector* list;
	float cutoff;
	string distFile;
	vector<int> rowPos;
	MothurOut* m;
};

//**********************************************************************************************************************

#endif


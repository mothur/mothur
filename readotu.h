#ifndef READOTU_H
#define READOTU_H
/*
 *  readotu.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */



#include "rabundvector.hpp"
#include "listvector.hpp"
#include "sharedlistvector.h"
#include "inputdata.h"
#include "globaldata.hpp"
#include "sabundvector.hpp"
#include "groupmap.h"


class ReadOTUFile {
	
public:
	ReadOTUFile(string);
	~ReadOTUFile();
	void read(GlobalData* globaldata);
private:
	//ifstream fileHandle;
	string philFile;
	InputData* input;
	InputData* inputSabund;
	InputData* inputRabund;
	InputData* inputList;
	ListVector* list;
	SharedListVector* SharedList;
	OrderVector* order;
	SAbundVector* sabund;
	RAbundVector* rabund;
	GlobalData* globaldata;
	// InputData* getInput()			{	return input;	}
};


#endif

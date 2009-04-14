#ifndef VENNCOMMAND_H
#define VENNCOMMAND_H
/*
 *  venncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "readmatrix.hpp"
#include "sharedlistvector.h"
#include "venn.h"
#include "validcalculator.h"
#include "sharedutilities.h"


class GlobalData;


class VennCommand : public Command {

public:
	VennCommand();
	~VennCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadMatrix* read;
	SharedUtil* util;
	InputData* input;
	SharedListVector* SharedList;
	SharedOrderVector* order;
	OrderVector* ordersingle;
	Venn* venn;
	string format;
	vector<Calculator*> vennCalculators;	
	ValidCalculators* validCalculator;
	int abund;

};



#endif
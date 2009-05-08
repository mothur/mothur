#ifndef DISTANCECOMMAND_H
#define DISTANCECOMMAND_H

/*
 *  distancecommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "globaldata.hpp"
#include "validcalculator.h"
#include "dist.h"
#include "sequencedb.h"


using namespace std;

class DistanceCommand : public Command {

public:
	DistanceCommand();	
	~DistanceCommand() {};
	int execute();	
	
private:
	GlobalData* globaldata;
	ValidCalculators* validCalculator;
	Dist* distCalculator;
	SequenceDB* seqDB;
	ofstream out;
	string outputFileName;
	string ends, distFile;
	int processors;
	float cutoff;
	
	int driver(Dist*, SequenceDB*, int, int, string, float);

};

#endif




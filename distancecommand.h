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
	ofstream out, phylipOut;
	ifstream in;
	string outputFileName;
	string countends;
	int processors;
	float cutoff;
	
	void appendFiles(string, string);
	int driver(Dist*, SequenceDB*, int, int, string, string, float);

};

#endif




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
#include "validcalculator.h"
#include "dist.h"
#include "sequencedb.h"


class DistanceCommand : public Command {

public:
	DistanceCommand(string);
	~DistanceCommand();
	int execute();	
	void help();
	
private:
	struct linePair {
		int start;
		int end;
	};
	
	Dist* distCalculator;
	SequenceDB alignDB;

	string countends, output, fastafile, calc, outputDir;
	int processors;
	float cutoff;
	map<int, int> processIDS;   //end line, processid
	vector<linePair*> lines;
	
	bool abort;
	vector<string>  Estimators; //holds estimators to be used
	
	//void appendFiles(string, string);
	void createProcesses(string);
	int driver(/*Dist*, SequenceDB, */int, int, string, float);
	int driverMPI(int, int, MPI_File&, float);
	
	int convertMatrix(string);
	int convertToLowerTriangle(string);

};

#endif




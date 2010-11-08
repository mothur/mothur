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
	DistanceCommand();
	~DistanceCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();	
	void help();
	
private:
	struct distlinePair {
		int start;
		int end;
		
	};
	
	Dist* distCalculator;
	SequenceDB alignDB;

	string countends, output, fastafile, calc, outputDir, oldfastafile, column, compress;

	int processors, numNewFasta;
	float cutoff;
	vector<int> processIDS;   //end line, processid
	vector<distlinePair> lines;
	
	bool abort;
	vector<string>  Estimators, outputNames; //holds estimators to be used
	map<string, vector<string> > outputTypes;
	
	//void m->appendFiles(string, string);
	void createProcesses(string);
	int driver(/*Dist*, SequenceDB, */int, int, string, float);
	int driver(int, int, string, string);
	
	#ifdef USE_MPI 
	int driverMPI(int, int, MPI_File&, float);
	int driverMPI(int, int, string, unsigned long int&);
	int driverMPI(int, int, string, unsigned long int&, string);
	#endif
	
	//int convertMatrix(string);
	bool sanityCheck();
	//int convertToLowerTriangle(string);

};

#endif




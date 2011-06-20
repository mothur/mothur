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
	~DistanceCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "dist.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";	}
	string getHelpString();	
	string getCitation() { return "Schloss PD (2010). The effects of alignment quality, distance calculation method, sequence filtering, and region on the analysis of 16S rRNA gene-based studies. PLoS Comput Biol 6: e1000844. \nhttp://www.mothur.org/wiki/Dist.seqs"; }
	string getDescription()		{ return "calculate the pairwaise distances between aligned sequences"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
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




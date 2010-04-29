#ifndef SEQSUMMARYCOMMAND_H
#define SEQSUMMARYCOMMAND_H

/*
 *  seqcoordcommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"

class SeqSummaryCommand : public Command {
public:
	SeqSummaryCommand(string);
	~SeqSummaryCommand();
	int execute();
	void help();
	
private:
	bool abort;
	string fastafile, outputDir;
	int processors;
	
	struct linePair {
		int start;
		int num;
		linePair(long int i, long int j) : start(i), num(j) {}
	};
	vector<linePair*> lines;
	vector<int> processIDS;
	
	int createProcessesCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string, string);
	int driverCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string, string, linePair*);	
	int setLines(string);

	#ifdef USE_MPI
	int MPICreateSummary(int, int, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, MPI_File&, MPI_File&, vector<long>&);	
	#endif


};

#endif

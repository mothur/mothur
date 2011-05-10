#ifndef FILTERSEQSCOMMAND_H
#define FILTERSEQSCOMMAND_H

/*
 *  filterseqscommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "filters.h"

class Sequence;
class FilterSeqsCommand : public Command {

public:
	FilterSeqsCommand(string);
	FilterSeqsCommand();
	~FilterSeqsCommand() {};
	
	vector<string> setParameters();
	string getCommandName()			{ return "filter.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Filter.seqs"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};

	vector<linePair*> lines;
	vector<int> processIDS;

	string vertical, filter, fasta, hard, outputDir, filterFileName;
	vector<string> fastafileNames;	
	int alignmentLength, processors;
	vector<int> bufferSizes;
	vector<string> outputNames;

	char trump;
	bool abort;
	float soft;
	int numSeqs;
	
	string createFilter();
	int filterSequences();
	int createProcessesCreateFilter(Filters&, string);
	int createProcessesRunFilter(string, string);
	int driverRunFilter(string, string, string, linePair*);
	int driverCreateFilter(Filters& F, string filename, linePair* line);
	#ifdef USE_MPI
	int driverMPIRun(int, int, MPI_File&, MPI_File&, vector<unsigned long int>&);
	int MPICreateFilter(int, int, Filters&, MPI_File&, vector<unsigned long int>&);	
	#endif
	
};

#endif

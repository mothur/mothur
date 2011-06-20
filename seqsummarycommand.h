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
	SeqSummaryCommand();
	~SeqSummaryCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.seqs"; }
	string getDescription()		{ return "summarize the quality of sequences in an unaligned or aligned fasta file"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:
	bool abort;
	string fastafile, outputDir, namefile;
	int processors;
	vector<string> outputNames;
	map<string, int> nameMap;
	
	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};

	vector<linePair*> lines;
	vector<int> processIDS;
	
	int createProcessesCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string, string);
	int driverCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string, string, linePair*);	

	#ifdef USE_MPI
	int MPICreateSummary(int, int, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, MPI_File&, MPI_File&, vector<unsigned long int>&);	
	#endif


};

#endif

#ifndef PAIRWISESEQSCOMMAND_H
#define PAIRWISESEQSCOMMAND_H

/*
 *  pairwiseseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "validcalculator.h"
#include "dist.h"
#include "sequencedb.h"

class PairwiseSeqsCommand : public Command {
	
public:
	PairwiseSeqsCommand(string);	
	PairwiseSeqsCommand();
	~PairwiseSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pairwise.seqs";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "Needleman SB, Wunsch CD (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol 48: 443-53. [ for needleman ]\nGotoh O (1982). An improved algorithm for matching biological sequences. J Mol Biol 162: 705-8. [ for gotoh ] \nhttp://www.mothur.org/wiki/Pairwise.seqs"; }
	string getDescription()		{ return "calculates pairwise distances from an unaligned fasta file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	struct distlinePair {
		int start;
		int end;
	};
	
	vector<int> processIDS;   //end line, processid
	vector<distlinePair> lines;
	
	Alignment* alignment;
	Dist* distCalculator;
	SequenceDB alignDB;
	
	void createProcesses(string);
	int driver(int, int, string, float);
	int driver(int, int, string, string);
	
	#ifdef USE_MPI 
	int driverMPI(int, int, MPI_File&, float);
	int driverMPI(int, int, string, unsigned long int&);
	int driverMPI(int, int, string, unsigned long int&, string);
	#endif
	
	string fastaFileName, align, calc, outputDir, output;
	float match, misMatch, gapOpen, gapExtend, cutoff;
	int processors;
	vector<string> fastaFileNames, Estimators;
	vector<string> outputNames;
	
	bool abort, countends, compress;
};

#endif


#ifndef ALIGNCOMMAND_H
#define ALIGNCOMMAND_H

/*
 *  aligncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "alignmentdb.h"

class AlignCommand : public Command {
	
public:
	AlignCommand(string);	
	~AlignCommand();
	int execute(); 
	void help();	
	
private:
	struct linePair {
		unsigned long int start;
		int numSeqs;
		linePair(unsigned long int i, int j) : start(i), numSeqs(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	bool MPIWroteAccnos;
	
	AlignmentDB* templateDB;
	Alignment* alignment;
	
	int driver(linePair*, string, string, string, string);
	int createProcesses(string, string, string, string);
	void appendAlignFiles(string, string); 
	void appendReportFiles(string, string);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&);
	#endif
	
	string candidateFileName, templateFileName, distanceFileName, search, align, outputDir;
	float match, misMatch, gapOpen, gapExtend, threshold;
	int processors, kmerSize;
	vector<string> candidateFileNames;
	
	bool abort, flip;
};

#endif

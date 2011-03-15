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
#include "mothurfiles.h"

class AlignCommand : public Command {
	
public:
	AlignCommand(string);	
	AlignCommand();
	~AlignCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute(); 
	void help();	
	
private:
	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	bool MPIWroteAccnos;
	map<string, vector<string> > outputTypes;
	
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
	vector<string> outputNames;
	
	bool abort, flip, calledHelp;
	CurrentFile* currentFiles;

};

#endif

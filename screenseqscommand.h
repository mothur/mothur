#ifndef SCREENSEQSCOMMAND_H
#define SCREENSEQSCOMMAND_H

/*
 *  screenseqscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/3/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */
#include "mothur.h"
#include "command.hpp"

class ScreenSeqsCommand : public Command {
	
public:
	ScreenSeqsCommand(string);
	~ScreenSeqsCommand();
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

	int screenNameGroupFile(set<string>);
	int screenGroupFile(set<string>);
	int screenAlignReport(set<string>);
	
	int driver(linePair*, string, string, string, string, set<string>&);
	int createProcesses(string, string, string, string, set<string>&);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&, set<string>&);
	#endif

	bool abort;
	string fastafile, namefile, groupfile, alignreport, outputDir;
	int startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength, processors;
	vector<string> outputNames;
};

#endif

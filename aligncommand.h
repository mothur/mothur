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

class AlignCommand : public Command {
	
public:
	AlignCommand(string);	
	~AlignCommand();
	int execute(); 
	void help();	
	
private:
	struct linePair {
		int start;
		int numSeqs;
		linePair(int i, int j) : start(i), numSeqs(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	
	Database* templateDB;
	Alignment* alignment;
	
	int driver(linePair*, string, string);
	void createProcesses(string, string);
	void appendAlignFiles(string, string); 
	void appendReportFiles(string, string);
	
	string candidateFileName, templateFileName, distanceFileName, search, align;
	float match, misMatch, gapOpen, gapExtend;
	int processors, kmerSize;
	
	bool abort;
};

#endif

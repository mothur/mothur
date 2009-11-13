#ifndef CLASSIFYSEQSCOMMAND_H
#define CLASSIFYSEQSCOMMAND_H

/*
 *  classifyseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/2/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "alignment.hpp"
#include "classify.h"

class ClassifySeqsCommand : public Command {
	
public:
	ClassifySeqsCommand(string);	
	~ClassifySeqsCommand();
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
	
	Classify* classify;
	
	string fastaFileName, templateFileName, distanceFileName, search, method, taxonomyFileName;
	int processors, kmerSize, numWanted, cutoff;
	float match, misMatch, gapOpen, gapExtend;
	bool abort;
	
	int driver(linePair*, string, string);
	void appendTaxFiles(string, string);
	void createProcesses(string, string); 
};

#endif


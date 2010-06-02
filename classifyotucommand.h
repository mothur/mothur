#ifndef CLASSIFYOTUSCOMMAND_H
#define CLASSIFYOTUSCOMMAND_H

/*
 *  classifyotucommand.h
 *  Mothur
 *
 *  Created by westcott on 6/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "listvector.hpp"
#include "inputdata.h"


class ClassifyOtuCommand : public Command {

public:
	ClassifyOtuCommand(string);
	~ClassifyOtuCommand();
	int execute();
	void help();

private:

	ListVector* list;
	InputData* input;
	string listfile, namefile, taxfile, label, outputDir;
	bool abort, allLines, probs;
	int cutoff;
	set<string> labels; //holds labels to be used
	vector<string> outputNames;
	map<string, string> nameMap;
	map<string, string> taxMap;

	int readNamesFile();
	int readTaxonomyFile();
	void removeConfidences(string&);
	int process(ListVector*);
	string findConsensusTaxonomy(int, ListVector*, int&); 	// returns the name of the "representative" taxonomy of given bin
	
												
};

#endif



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
	ClassifyOtuCommand();
	~ClassifyOtuCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "classify.otu";		}
	string getCommandCategory()		{ return "Phylotype Analysis";	}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "Schloss PD, Westcott SL (2011). Assessing and improving methods used in OTU-based approaches for 16S rRNA gene sequence analysis. Appl Environ Microbiol 77:3219.\nhttp://www.mothur.org/wiki/Classify.otu"; }
	string getDescription()		{ return "find the concensus taxonomy for each OTU"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	

private:

	ListVector* list;
	InputData* input;
	string listfile, namefile, taxfile, label, outputDir, refTaxonomy, groupfile, basis;
	bool abort, allLines, probs;
	int cutoff;
	set<string> labels; //holds labels to be used
	vector<string> outputNames;
	map<string, string> nameMap;
	map<string, string> taxMap;

	int process(ListVector*);
	string addUnclassifieds(string, int);
	vector<string> findConsensusTaxonomy(int, ListVector*, int&, string&); 	// returns the name of the "representative" taxonomy of given bin
	
												
};

#endif



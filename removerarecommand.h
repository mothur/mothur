#ifndef REMOVERARECOMMAND_H
#define REMOVERARECOMMAND_H

/*
 *  removerarecommand.h
 *  mothur
 *
 *  Created by westcott on 1/21/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "listvector.hpp"

class RemoveRareCommand : public Command {
	
public:
	
	RemoveRareCommand(string);	
	RemoveRareCommand();
	~RemoveRareCommand(){}
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	
	
private:
	string sabundfile, rabundfile, sharedfile, groupfile, listfile, outputDir, groups, label;
	int nseqs, allLines;
	bool abort, byGroup;
	vector<string> outputNames, Groups;
	set<string> labels;
	map<string, vector<string> > outputTypes;
	
	int processSabund();
	int processRabund();
	int processList();
	int processShared();
	
};

#endif





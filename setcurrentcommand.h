#ifndef SETCURRENTCOMMAND_H
#define SETCURRENTCOMMAND_H

/*
 *  setcurrentcommand.h
 *  Mothur
 *
 *  Created by westcott on 3/16/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"

class SetCurrentCommand : public Command {
	
public:
	SetCurrentCommand(string);
	SetCurrentCommand();
	~SetCurrentCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	bool abort;
	
	string clearTypes;
	vector<string> types;
	
	string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile;
	string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile;

	
};

#endif



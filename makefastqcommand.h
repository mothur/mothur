#ifndef MAKEFASTQCOMMAND_H
#define MAKEFASTQCOMMAND_H

/*
 *  makefastqcommand.h
 *  mothur
 *
 *  Created by westcott on 2/14/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"

class MakeFastQCommand : public Command {
	
public:
	
	MakeFastQCommand(string);	
	MakeFastQCommand();
	~MakeFastQCommand(){}
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	
	
private:
	
	string fastafile, qualfile, outputDir;
	bool abort;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	
	string convertQual(vector<int>);
	
};

#endif



#ifndef SHAREDCOMMAND_H
#define SHAREDCOMMAND_H
/*
 *  sharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "sharedlistvector.h"
#include "inputdata.h"
#include "readotu.h"

/* The shared() command:
	The shared command can only be executed after a successful read.shared command.  
	The shared command parses a .list file and separates it into groups.  
	It outputs a .shared file containing the OTU information for each group.  
	There are no shared command parameters.  The shared command should be in the following format: shared(). */


class GlobalData;

class SharedCommand : public Command {
	
public:
	SharedCommand(string);	
	SharedCommand();
	~SharedCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();	
	void help() {}
	
private:
	void printSharedData(vector<SharedRAbundVector*>);
	int createMisMatchFile();
	int readOrderFile();
	bool isValidGroup(string, vector<string>);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
	map<string, vector<string> > outputTypes;
	
	GlobalData* globaldata;
	ReadOTUFile* read;
	SharedListVector* SharedList;
	InputData* input;
	GroupMap* groupMap;
	vector<string> groups, outputNames, order;
	ofstream out;
	string filename, fileroot, outputDir;
	bool firsttime, pickedGroups;
	map<string, ofstream*> filehandles;
	map<string, ofstream*>::iterator it3;

};

#endif

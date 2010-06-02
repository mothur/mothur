#ifndef GETLISTCOUNTCOMMAND_H
#define GETLISTCOUNTCOMMAND_H
/*
 *  getlistcountcommand.h
 *  Mothur
 *
 *  Created by westcott on 10/12/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "listvector.hpp"
#include "readotu.h"

class GlobalData;

/**********************************************************/

class GetListCountCommand : public Command {
	
public:
	GetListCountCommand(string);	
	~GetListCountCommand();
	int execute();
	void help();	
	
private:
	GlobalData* globaldata;
	ListVector* list;
	ReadOTUFile* read;
	InputData* input;
	
	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label, listfile, outputDir, sort;
	ofstream out;
	vector<string> outputNames;
	
	void process(ListVector*);
};
/**********************************************************/

#endif

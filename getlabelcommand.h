#ifndef GETLABELCOMMAND_H
#define GETLABELCOMMAND_H

/*
 *  getlabelcommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "ordervector.hpp"
#include "inputdata.h"


class GetlabelCommand : public Command {
public:
	GetlabelCommand(string);
	GetlabelCommand();
	~GetlabelCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.label";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
    string getOutputPattern(string) { return ""; }	
	string getCitation() { return "http://www.mothur.org/wiki/Get.label"; }
	string getDescription()		{ return "outputs labels"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	string inputfile, listfile, rabundfile, sabundfile, format;
	bool abort;
	vector<string> outputNames;
};

#endif

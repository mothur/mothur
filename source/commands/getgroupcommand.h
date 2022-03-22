#ifndef GETGROUPCOMMAND_H
#define GETGROUPCOMMAND_H

/*
 *  getgroupcommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"

class GetgroupCommand : public Command {
public:
	GetgroupCommand(string);
	~GetgroupCommand() = default;
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.group";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
    string getOutputPattern(string)	{ return "";  }
	string getCitation() { return "http://www.mothur.org/wiki/Get.group"; }
	string getDescription()		{ return "outputs group names"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	
	string outputFile, sharedfile;
	vector<string> outputNames;
	ofstream out;
	ifstream in;
	bool abort;

};

#endif

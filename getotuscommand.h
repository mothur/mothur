#ifndef GETOTUSCOMMAND_H
#define GETOTUSCOMMAND_H

/*
 *  getotuscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "command.hpp"
#include "groupmap.h"
#include "listvector.hpp"

class GetOtusCommand : public Command {
	
public:
	
	GetOtusCommand(string);	
	GetOtusCommand();
	~GetOtusCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.otus";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Get.otus"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	string accnosfile, groupfile, listfile, outputDir, groups, label;
	bool abort;
	vector<string> outputNames, Groups;
	GroupMap* groupMap;
	
	void readAccnos();
	int readListGroup();
	int processList(ListVector*&, GroupMap*&, ofstream&, ofstream&, bool&);
	
};

#endif




#ifndef READTREECOMMAND_H
#define READTREECOMMAND_H

/*
 *  readtreecommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "readtree.h"
#include "treemap.h"


class ReadTreeCommand : public Command {
public:
	ReadTreeCommand(string);
	ReadTreeCommand() { abort = true; calledHelp = true; }
	~ReadTreeCommand() {}
	
	vector<string> setParameters() {  return outputNames; } //dummy doesn't really do anything
	string getCommandName()			{ return "read.tree";	}
	string getCommandCategory()		{ return "Hidden";	}
	string getHelpString() { return "This command is no longer available. You can provide your files directly to the downstream commands like unifrac.unweighted."; }	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	ReadTree* read;
	TreeMap* treeMap;
	string filename, treefile, groupfile, namefile;
	bool abort;
	map<string, string> nameMap;
	vector<string> outputNames;
	
	int readNamesFile();
	int numUniquesInName;

};


#endif

#ifndef DEUNIQUETREECOMMAND_H
#define DEUNIQUETREECOMMAND_H

/*
 *  deuniquetreecommand.h
 *  Mothur
 *
 *  Created by westcott on 5/27/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sharedutilities.h"
#include "readtree.h"

class DeuniqueTreeCommand : public Command {
	
public:
	DeuniqueTreeCommand(string);	
	DeuniqueTreeCommand();
	~DeuniqueTreeCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "deunique.tree";		}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Deunique.tree"; }
	string getDescription()		{ return "add the redundant sequence names back into a tree of unique sequences"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
	
private:
	int numUniquesInName;
	
	bool abort;
	string outputDir, treefile, namefile;
	vector<string> outputNames;
	map<string, string> nameMap;
	int readNamesFile();
};

#endif

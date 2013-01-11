#ifndef PHYLOTYPECOMMAND_H
#define PHYLOTYPECOMMAND_H


/*
 *  phylotypecommand.h
 *  Mothur
 *
 *  Created by westcott on 11/20/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "taxonomyequalizer.h"
#include "command.hpp"

/*************************************************************************/

class PhylotypeCommand : public Command {
	
public:
	PhylotypeCommand(string);	
	PhylotypeCommand();
	~PhylotypeCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "phylotype";		}
	string getCommandCategory()		{ return "Clustering";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Phylotype"; }
	string getDescription()		{ return "cluster your sequences into OTUs based on their classifications"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:
	bool abort, allLines;
	string taxonomyFileName, label, outputDir, namefile;
	set<string> labels; //holds labels to be used
	int cutoff;
	map<string, string> namemap;
	vector<string> outputNames;
	
	map<int, int> currentNodes;
	map<int, int> parentNodes;
	map<int, int>::iterator itCurrent;
	
	int readNamesFile();

};


/*************************************************************************/


#endif




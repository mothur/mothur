#ifndef OTUHIERARCHYCOMMAND_H
#define OTUHIERARCHYCOMMAND_H
/*
 *  otuhierarchycommand.h
 *  Mothur
 *
 *  Created by westcott on 1/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "listvector.hpp"

//**********************************************************************************************************************

class OtuHierarchyCommand : public Command {

public:
	OtuHierarchyCommand(string);
	OtuHierarchyCommand();
	~OtuHierarchyCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "otu.hierarchy";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Otu.hierarchy"; }
	string getDescription()		{ return "relates OTUs at different distances"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	bool abort;
	set<string> labels; //holds labels to be used
	string label, listFile, outputDir, output;
	vector<string> outputNames;
	
	vector<ListVector> getListVectors();
		
};

//**********************************************************************************************************************

#endif



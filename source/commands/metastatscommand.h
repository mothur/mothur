#ifndef METASTATSCOMMAND_H
#define METASTATSCOMMAND_H

/*
 *  metastatscommand.h
 *  Mothur
 *
 *  Created by westcott on 9/16/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "sharedrabundvectors.hpp"
#include "mothurmetastats.h"
#include "designmap.h"


/**************************************************************************************************/

class MetaStatsCommand : public Command {

public:
	MetaStatsCommand(string);
	MetaStatsCommand();
	~MetaStatsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "metastats";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "White JR, Nagarajan N, Pop M (2009). Statistical methods for detecting differentially abundant features in clinical metagenomic samples. PLoS Comput Biol 5: e1000352. \nhttp://www.mothur.org/wiki/Metastats"; }
	string getDescription()		{ return "detects differentially abundant features in clinical metagenomic samples"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    DesignMap* designMap;
	SharedRAbundVectors* lookup;
		
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, inputDir, designfile, sets, sharedfile;
	vector<string> Groups, outputNames, Sets;
	vector< vector<string> > namesOfGroupCombos;
	int iters, processors;
	float threshold;
	
	int process(SharedRAbundVectors*&);
};

/**************************************************************************************************/

#endif


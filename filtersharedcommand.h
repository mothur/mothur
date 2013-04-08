//
//  filtersharedcommand.h
//  Mothur
//
//  Created by Sarah Westcott on 1/4/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_filtersharedcommand_h
#define Mothur_filtersharedcommand_h

#include "command.hpp"
#include "sharedrabundvector.h"
#include "inputdata.h"


class FilterSharedCommand : public Command {
    
public:
	FilterSharedCommand(string);
	FilterSharedCommand();
	~FilterSharedCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "filter.shared";	}
	string getCommandCategory()		{ return "OTU-Based Approaches";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Filter.shared"; }
	string getDescription()		{ return "remove OTUs based on various criteria"; }
    
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:	
	bool abort, pickedGroups, allLines, makeRare, keepties;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, sharedfile;
	vector<string> Groups, outputNames;
	int minAbund, minTotal, minSamples;
    float minPercent, minPercentSamples, rarePercent;
    
    int processShared(vector<SharedRAbundVector*>&);
	
};



#endif

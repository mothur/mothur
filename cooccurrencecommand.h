#ifndef COOCCURRENCECOMMAND_H
#define COOCCURRENCECOMMAND_H

/*
 *  COOCCURRENCE.h
 *  Mothur
 *
 *  Created by westcott on 11/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "trialswap2.h"
#include "inputdata.h"
#include "sharedrabundvector.h"


class CooccurrenceCommand : public Command {
	
public:
	
	CooccurrenceCommand(string);	
	CooccurrenceCommand();
	~CooccurrenceCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "Cooccurrence";			}
	string getCommandCategory()		{ return "Hypothesis Testing";	}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "Ulrich W & Gotelli NJ (2010).  Null model analysis of species associations using abundance data.  Ecology  91:3384.\nhttp://www.mothur.org/wiki/Cooccurrence"; }
	string getDescription()		{ return "calculates four metrics and tests their significance to assess whether presence-absence patterns are different than what one would expect by chance."; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
    string metric, matrix, outputDir;
    string label, sharedfile, groups;
    bool abort, allLines;
    set<string> labels;
    vector<string> outputNames, Groups;
    int runs;
    
    int getCooccurrence(vector<SharedRAbundVector*>&, ofstream&);
	
};

#endif



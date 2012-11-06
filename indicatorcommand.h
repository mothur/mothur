#ifndef INDICATORCOMMAND_H
#define INDICATORCOMMAND_H

/*
 *  indicatorcommand.h
 *  Mothur
 *
 *  Created by westcott on 11/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "readtree.h"
#include "counttable.h"
#include "sharedrabundvector.h"
#include "sharedrabundfloatvector.h"
#include "inputdata.h"

class IndicatorCommand : public Command {
public:
	IndicatorCommand(string);
	IndicatorCommand();
	~IndicatorCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "indicator";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "Dufrene M, Legendre P (1997). Species assemblages and indicator species: The need for a flexible asymmetrical approach. Ecol Monogr 67: 345-66.\n McCune B, Grace JB, Urban DL (2002). Analysis of ecological communities. MjM Software Design: Gleneden Beach, OR. \nLegendre P, Legendre L (1998). Numerical Ecology. Elsevier: New York. \nhttp://www.mothur.org/wiki/Indicator"; }
	string getDescription()		{ return "calculate the indicator value for each OTU"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	ReadTree* read;
	CountTable* ct;
	GroupMap* designMap;
	string treefile, sharedfile, relabundfile, groups, label, inputFileName, outputDir, designfile;
	bool abort;
	int iters, processors;
	vector<string> outputNames, Groups;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	
	int getShared();
	int getSharedFloat();
	int GetIndicatorSpecies(Tree*&);
	int GetIndicatorSpecies();
	set<string> getDescendantList(Tree*&, int, map<int, set<string> >, map<int, set<int> >&);
	vector<float> getValues(vector< vector<SharedRAbundVector*> >&, vector<string>&, map< vector<int>, vector<int> >);
	vector<float> getValues(vector< vector<SharedRAbundFloatVector*> >&, vector<string>&, map< vector<int>, vector<int> >);
    
	map<int, float> getDistToRoot(Tree*&);
	map< vector<int>, vector<int> > randomizeGroupings(vector< vector<SharedRAbundVector*> >&, int);
	map< vector<int>, vector<int> > randomizeGroupings(vector< vector<SharedRAbundFloatVector*> >&, int);
    
	vector<float> driver(vector< vector<SharedRAbundFloatVector*> >&, map< vector<int>, vector<int> >, int, vector<float>, int);
	vector<float> driver(vector< vector<SharedRAbundVector*> >&, map< vector<int>, vector<int> >, int, vector<float>, int);
    
	vector<float> getPValues(vector< vector<SharedRAbundFloatVector*> >&, map< vector<int>, vector<int> >, int, vector<float>);
	vector<float> getPValues(vector< vector<SharedRAbundVector*> >&, map< vector<int>, vector<int> >, int, vector<float>);

	
};


#endif


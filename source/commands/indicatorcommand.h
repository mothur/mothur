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


#include "inputdata.h"
#include "designmap.h"

//**********************************************************************************************************************
struct sharedIndexes {
    int groupIndex;
    int otuIndex;
    
    sharedIndexes() : groupIndex(0), otuIndex(0) {}
    sharedIndexes(int g, int o) : groupIndex(g), otuIndex(o) {}
    
    bool operator<(const sharedIndexes& rhs) const {
        return rhs.groupIndex < this->groupIndex || (rhs.groupIndex == this->groupIndex && rhs.otuIndex < this->otuIndex);
    }
    
    bool operator>(const sharedIndexes& rhs) const {
        return rhs.groupIndex > this->groupIndex || (rhs.groupIndex == this->groupIndex && rhs.otuIndex > this->otuIndex);
    }
    
    bool operator=(const sharedIndexes& rhs) const {
        return ((rhs.groupIndex == this->groupIndex) && (rhs.otuIndex == this->otuIndex));
    }
};
//**********************************************************************************************************************

class IndicatorCommand : public Command {
public:
	IndicatorCommand(string);
	~IndicatorCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "indicator";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Dufrene M, Legendre P (1997). Species assemblages and indicator species: The need for a flexible asymmetrical approach. Ecol Monogr 67: 345-66.\n McCune B, Grace JB, Urban DL (2002). Analysis of ecological communities. MjM Software Design: Gleneden Beach, OR. \nLegendre P, Legendre L (1998). Numerical Ecology. Elsevier: New York. \nhttp://www.mothur.org/wiki/Indicator"; }
	string getDescription()		{ return "calculate the indicator value for each OTU"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	DesignMap* designMap;
	string treefile, sharedfile, relabundfile, groups, label, inputFileName,  designfile;
    bool abort;
    int iters, processors;
    vector<string> outputNames, Groups, namesSeqs;
    	
	set<string> getDescendantList(Tree*&, int, map<int, set<string> >, map<int, set<int> >&);
	map<int, float> getDistToRoot(Tree*&);
	vector< map<sharedIndexes, sharedIndexes> > randomizeGroupings(vector<int>, int);
    
    SharedRAbundVectors* getShared();
    SharedRAbundFloatVectors* getSharedFloat();
    int GetIndicatorSpecies(Tree*&);
    int GetIndicatorSpecies();
	vector<float> getPValues(vector< vector<SharedRAbundFloatVector*> >&, vector< vector<string> >&, int, vector<float>);
	vector<float> getPValues(vector< vector<SharedRAbundVector*> >&, vector< vector<string> >&, int, vector<float>);
};
/**************************************************************************************************/

#endif


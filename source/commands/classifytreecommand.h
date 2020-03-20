#ifndef Mothur_classifytreecommand_h
#define Mothur_classifytreecommand_h

//
//  classifytreecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 2/20/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "command.hpp"
#include "readtree.h"
#include "treemap.h"
#include "counttable.h"

class ClassifyTreeCommand : public Command {
public:
	ClassifyTreeCommand(string);
	~ClassifyTreeCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "classify.tree";				}
	string getCommandCategory()		{ return "Phylotype Analysis";          }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Classify.tree"; }
	string getDescription()		{ return "Find the consensus taxonomy for the descendant of each tree node"; }
    
	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string treefile, taxonomyfile, groupfile, namefile, countfile, outputDir, output;
	bool abort;
	vector<string> outputNames;
    int numUniquesInName, cutoff;
    map<string, string> nameMap;
    map<string, int> nameCount;
    map<string, string> taxMap;
    CountTable* ct;
	
	int getClassifications(Tree*&);
	map<string, set<string> > getDescendantList(Tree*&, int, map<int, map<string, set<string> > >);
    string getTaxonomy(set<string>, int&);
	
};



#endif

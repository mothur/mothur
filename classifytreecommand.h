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

class ClassifyTreeCommand : public Command {
public:
	ClassifyTreeCommand(string);
	ClassifyTreeCommand();
	~ClassifyTreeCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "classify.tree";				}
	string getCommandCategory()		{ return "Phylotype Analysis";          }
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Classify.tree"; }
	string getDescription()		{ return "Find the consensus taxonomy for the descendant of each tree node"; }
    
	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	ReadTree* read;
    TreeMap* tmap;
	string treefile, taxonomyfile, groupfile, namefile, outputDir;
	bool abort;
	vector<string> outputNames;
    int numUniquesInName, cutoff;
    map<string, string> nameMap;
    map<string, int> nameCount;
    map<string, string> taxMap;
	
	int getClassifications(Tree*&);
	map<string, set<string> > getDescendantList(Tree*&, int, map<int, map<string, set<string> > >);
    string getTaxonomy(set<string>, int&);
    int readNamesFile(); 
    int readTaxonomyFile();
	
};



#endif

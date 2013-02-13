//
//  mergetaxsummarycommand.h
//  Mothur
//
//  Created by Sarah Westcott on 2/13/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_mergetaxsummarycommand_h
#define Mothur_mergetaxsummarycommand_h

#include "mothur.h"
#include "command.hpp"
#include "phylosummary.h"

class MergeTaxSummaryCommand : public Command {
public:
	MergeTaxSummaryCommand(string);
	MergeTaxSummaryCommand();
	~MergeTaxSummaryCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "merge.taxsummary";	}
	string getCommandCategory()		{ return "Phylotype Analysis";		}
	string getHelpString();	
    string getOutputPattern(string){ return "";  }	
	string getCitation() { return "http://www.mothur.org/wiki/Merge.taxsummary"; }
	string getDescription()		{ return "merges tax summary files creating one file"; }
    
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	vector<string> fileNames, outputNames;
	string outputFileName;
	int numInputFiles;
	bool abort;
       
    int addTaxToTree(vector<rawTaxNode>&, int, int, string, int, map<string, int>);
    int assignRank(int index, vector<rawTaxNode>& tree);
    int print(ofstream& out, vector<rawTaxNode>& tree, set<string> groups);
    int print(int, ofstream& out, vector<rawTaxNode>& tree, set<string> groups);
};


#endif

#ifndef SUMMARYTAXCOMMAND_H
#define SUMMARYTAXCOMMAND_H

/*
 *  summarytaxcommand.h
 *  Mothur
 *
 *  Created by westcott on 9/23/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "counttable.h"

/**************************************************************************************************/

class SummaryTaxCommand : public Command {
	public:
		SummaryTaxCommand(string);
		~SummaryTaxCommand(){}
		
		vector<string> setParameters();
		string getCommandName()			{ return "summary.tax";			}
		string getCommandCategory()		{ return "Phylotype Analysis";		}
		
	string getHelpString();	
    string getOutputPattern(string);	
		string getCitation() { return "http://www.mothur.org/wiki/Summary.tax"; }
		string getDescription()		{ return "summarize the taxonomies of a set of sequences"; }
		
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
		
	private:
		bool abort, relabund;
		string taxfile, outputDir, namefile, groupfile, countfile, output;
        int printlevel, threshold;
		vector<string> outputNames;
		map<string, int> nameMap;
    
        string processTaxMap(string);
        int findMaxLevel(string);
};

/**************************************************************************************************/


#endif

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
		SummaryTaxCommand();
		~SummaryTaxCommand(){}
		
		vector<string> setParameters();
		string getCommandName()			{ return "summary.tax";			}
		string getCommandCategory()		{ return "Phylotype Analysis";		}
		string getOutputFileNameTag(string, string);
	string getHelpString();	
		string getCitation() { return "http://www.mothur.org/wiki/Summary.tax"; }
		string getDescription()		{ return "summarize the taxonomies of a set of sequences"; }
		
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
		
	private:
		bool abort;
		string taxfile, outputDir, namefile, groupfile, refTaxonomy, countfile;
		vector<string> outputNames;
		map<string, int> nameMap;
};

/**************************************************************************************************/


#endif

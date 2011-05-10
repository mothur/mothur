#ifndef GETLINEAGECOMMAND_H
#define GETLINEAGECOMMAND_H

/*
 *  getlineagecommand.h
 *  Mothur
 *
 *  Created by westcott on 9/24/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"

class GetLineageCommand : public Command {
	
	public:
	
		GetLineageCommand(string);
		GetLineageCommand();
		~GetLineageCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "get.lineage";				}
		string getCommandCategory()		{ return "Phylotype Analysis";		}
		string getHelpString();	
		string getCitation() { return "http://www.mothur.org/wiki/Get.lineage"; }
	
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
		set<string> names;
		vector<string> outputNames;
		string fastafile, namefile, groupfile, alignfile, listfile, taxfile, outputDir, taxons;
		bool abort, dups;
		
		int readFasta();
		int readName();
		int readGroup();
		int readAlign();
		int readList();
		int readTax();	
		string removeConfidences(string);
		vector< map<string, float> > getTaxons(string);
};

#endif


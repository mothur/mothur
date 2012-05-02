#ifndef REMOVELINEAGECOMMAND_H
#define REMOVELINEAGECOMMAND_H

/*
 *  removelineagecommand.h
 *  Mothur
 *
 *  Created by westcott on 9/24/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"

class RemoveLineageCommand : public Command {
	
	public:
	
		RemoveLineageCommand(string);
		RemoveLineageCommand();
		~RemoveLineageCommand(){};
	
		vector<string> setParameters();
		string getCommandName()			{ return "remove.lineage";			}
		string getCommandCategory()		{ return "Phylotype Analysis";		}
		string getHelpString();	
		string getCitation() { return "http://www.mothur.org/wiki/Remove.lineage"; }
		string getDescription()		{ return "removes sequences from a list, fasta, name, group, alignreport or taxonomy file from a given taxonomy or set of taxonomies"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	private:
		set<string> names;
		vector<string> outputNames, listOfTaxons;
		string fastafile, namefile, groupfile, alignfile, listfile, taxfile, outputDir, taxons;
		bool abort, dups;
		
		int readFasta();
		int readName();
		int readGroup();
		int readAlign();
		int readList();
		int readTax();	
		vector< map<string, float> > getTaxons(string);
};

#endif


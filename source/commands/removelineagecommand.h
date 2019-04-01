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
#include "sharedrabundvectors.hpp"
#include "listvector.hpp"
#include "removeseqscommand.h"
#include "removeotulabelscommand.h"

class RemoveLineageCommand : public Command {
	
	public:
	
		RemoveLineageCommand(string);
		RemoveLineageCommand();
		~RemoveLineageCommand(){};
	
		vector<string> setParameters();
		string getCommandName()			{ return "remove.lineage";			}
		string getCommandCategory()		{ return "Phylotype Analysis";		}
		
	string getHelpString();	
    string getOutputPattern(string);	
		string getCitation() { return "http://www.mothur.org/wiki/Remove.lineage"; }
		string getDescription()		{ return "removes sequences from a list, fasta, name, group, alignreport or taxonomy file from a given taxonomy or set of taxonomies"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	private:
		vector<string> outputNames, listOfTaxons;
		string fastafile, namefile, groupfile, alignfile, listfile, countfile, taxfile, outputDir, taxons, sharedfile, constaxonomy, label;
		bool abort, dups;
		
        string readTax();
        string readConsTax();
        int runRemoveOTUs(string);
        int runRemoveSeqs(string);
		vector< map<string, float> > getTaxons(string);
};

#endif


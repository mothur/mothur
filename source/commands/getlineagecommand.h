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
#include "sharedrabundvectors.hpp"
#include "listvector.hpp"
#include "getseqscommand.h"
#include "getotulabelscommand.h"

class GetLineageCommand : public Command {
	
	public:
	
		GetLineageCommand(string);
		GetLineageCommand();
		~GetLineageCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "get.lineage";				}
		string getCommandCategory()		{ return "Phylotype Analysis";		}
		
	string getHelpString();	
    string getOutputPattern(string);	
		string getCitation() { return "http://www.mothur.org/wiki/Get.lineage"; }
		string getDescription()		{ return "gets sequences from a list, fasta, name, group, alignreport or taxonomy file from a given taxonomy or set of taxonomies"; }

	
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
		vector<string> outputNames, listOfTaxons;
		string fastafile, namefile, groupfile, alignfile, countfile, listfile, taxfile, outputDir, taxons, sharedfile, constaxonomy, label, accnosFileName;
		bool abort, dups;
		
		int readTax();
        int readConsTax();
    
        int runGetOTUs();
        int runGetSeqs();
    
    
		
};

#endif


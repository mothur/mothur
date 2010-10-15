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
		~GetLineageCommand(){};
		vector<string> getRequiredParameters();
		vector<string> getValidParameters();
		vector<string> getRequiredFiles();
		map<string, vector<string> > getOutputFiles() { return outputTypes; }
		int execute();
		void help();	
		
	private:
		set<string> names;
		vector<string> outputNames;
		string fastafile, namefile, groupfile, alignfile, listfile, taxfile, outputDir, taxons;
		bool abort, dups;
		map<string, vector<string> > outputTypes;
		
		int readFasta();
		int readName();
		int readGroup();
		int readAlign();
		int readList();
		int readTax();	
		string removeConfidences(string);
};

#endif


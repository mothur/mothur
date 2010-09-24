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
		~RemoveLineageCommand(){};
		int execute();
		void help();	
		
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
};

#endif


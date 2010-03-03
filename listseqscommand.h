#ifndef LISTSEQSCOMMAND_H
#define LISTSEQSCOMMAND_H

/*
 *  listseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"

class ListSeqsCommand : public Command {
	
	public:
	
		ListSeqsCommand(string);	
		~ListSeqsCommand(){};
		int execute();
		void help();	
		
	private:
		vector<string> names;
		string fastafile, namefile, groupfile, alignfile, inputFileName, outputDir, listfile;
		bool abort;
		
		int readFasta();
		int readName();
		int readGroup();
		int readAlign();
		int readList();
		
};

#endif


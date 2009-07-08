#ifndef GETSEQSCOMMAND_H
#define GETSEQSCOMMAND_H

/*
 *  getseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "command.hpp"

class GetSeqsCommand : public Command {
	
	public:
	
		GetSeqsCommand(string);	
		~GetSeqsCommand(){};
		int execute();
		void help();	
		
	private:
		set<string> names;
		string accnosfile, fastafile, namefile, groupfile, alignfile;
		bool abort;
		
		void readFasta();
		void readName();
		void readGroup();
		void readAlign();
		void readAccnos();
		
};

#endif


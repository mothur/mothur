#ifndef REMOVESEQSCOMMAND_H
#define REMOVESEQSCOMMAND_H

/*
 *  removeseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"

class RemoveSeqsCommand : public Command {
	
	public:
	
		RemoveSeqsCommand(string);	
		RemoveSeqsCommand();
		~RemoveSeqsCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "remove.seqs";				}
		string getCommandCategory()		{ return "Sequence Processing";		}
		string getHelpString();	
	
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
		set<string> names;
		string accnosfile, fastafile, namefile, groupfile, alignfile, listfile, taxfile, qualfile, outputDir;
		bool abort, dups;
		vector<string> outputNames;
		
		int readFasta();
		int readName();
		int readGroup();
		int readAlign();
		void readAccnos();
		int readList();
		int readTax();
		int readQual();
		
};

#endif


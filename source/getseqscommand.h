#ifndef GETSEQSCOMMAND_H
#define GETSEQSCOMMAND_H

/*
 *  getseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"

class GetSeqsCommand : public Command {
	
	public:
	
		GetSeqsCommand(string);
		GetSeqsCommand();
		~GetSeqsCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "get.seqs";				}
		string getCommandCategory()		{ return "Sequence Processing";		}
		string getHelpString();	
		string getCitation() { return "http://www.mothur.org/wiki/Get.seqs"; }
		string getDescription()		{ return "gets sequences from a list, fasta, name, group, alignreport, quality or taxonomy file"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
		set<string> names;
		vector<string> outputNames;
		string accnosfile, accnosfile2, fastafile, namefile, groupfile, alignfile, listfile, taxfile, qualfile, outputDir;
		bool abort, dups;
		
		int readFasta();
		int readName();
		int readGroup();
		int readAlign();
		int readAccnos();
		int readList();
		int readTax();
		int readQual();
		int compareAccnos();
		
};

#endif


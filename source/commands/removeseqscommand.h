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
		~RemoveSeqsCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "remove.seqs";				}
		string getCommandCategory()		{ return "Sequence Processing";		}
		
	string getHelpString();	
    string getOutputPattern(string);	
		string getCitation() { return "http://www.mothur.org/wiki/Remove.seqs"; }
		string getDescription()		{ return "removes sequences from a list, fasta, name, group, alignreport, contigsreport, quality or taxonomy file"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	private:
		set<string> names;
		string accnosfile, fastafile, fastqfile, namefile, groupfile, countfile, alignfile, listfile, taxfile, qualfile,  format, contigsreportfile;
		bool abort, dups;
		vector<string> outputNames;
        map<string, string> uniqueMap;
		
		void readFasta();
        void readFastq();
		void readName();
		void readGroup();
        void readCount();
		void readAlign();
        void readContigs();
		void readList();
		void readTax();
		void readQual();
		
};

#endif


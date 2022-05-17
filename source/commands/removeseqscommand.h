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
        RemoveSeqsCommand(string, pair<string, string> dupsFile, string dupsFileType);
        RemoveSeqsCommand(unordered_set<string>, pair<string, string> dupsFile, string dupsFileType);
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
        unordered_set<string> names;
        vector<string> fastafiles, namefiles, groupfiles, countfiles, alignfiles, listfiles, taxfiles, fastqfiles, contigsreportfiles, qualityfiles, outputNames;
		string accnosfile, format;
		bool abort, dups;
        map<string, string> uniqueMap;
		
		void readFasta(string);
        void readFastq(string);
        void readGZFastq(string);
        void readName(string); //inputNameFile, mothur generates output name
        void readName(string, string); //inputNameFile, outputName (internal use)
		void readGroup(string);
        void readCount(string); //inputCountFile, mothur generates output name
        void readCount(string, string); //inputCountFile, outputName (internal use)
		void readAlign(string);
        void readContigs(string);
		void readList(string);
		void readTax(string);
		void readQual(string);
		
};

#endif


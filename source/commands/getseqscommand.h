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
    string getOutputPattern(string);	
		string getCitation() { return "http://www.mothur.org/wiki/Get.seqs"; }
		string getDescription()		{ return "gets sequences from a list, fasta, count, name, group, alignreport, quality, fastq, contigsreport or taxonomy file"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
		set<string> names;
		vector<string> outputNames;
		string accnosfile, accnosfile2, fastafile, fastqfile, namefile, countfile, groupfile, alignfile, listfile, taxfile, qualfile, outputDir, format, contigsreportfile;
		bool abort, dups;
        map<string, string> uniqueMap;
        //for debug
        map<string, set<string> > sanity; //maps file type to names chosen for file. something like "fasta" -> vector<string>. If running in debug mode this is filled and we check to make sure all the files have the same names. If they don't we output the differences for the user.
		
		int readFasta();
        int readFastq();
		int readName();
		int readGroup();
        int readCount();
		int readAlign();
		int readList();
		int readTax();
		int readQual();
        int readContigsReport();
		int compareAccnos();
        int runSanityCheck();
        int createMisMatchFile(ofstream&, string, string, set<string>, set<string>);

		
};

#endif


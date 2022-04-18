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
#include "sequencedb.h"

//**********************************************************************************************************************


class GetSeqsCommand : public Command {
	
	public:
	
		GetSeqsCommand(string);
        GetSeqsCommand(unordered_set<string>, pair<string,string> fasta, pair<string,vector<string> > list, pair<string,string> dupsFile, string dupsFileType);
        GetSeqsCommand(unordered_set<string>, pair<string,string> fasta, pair<string,string> list, pair<string,string> dupsFile, string dupsFileType);
        GetSeqsCommand(map<string, vector<int> >, string fastafile, vector<string> outputFiles, vector<string> groups);
		~GetSeqsCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "get.seqs";				}
		string getCommandCategory()		{ return "Sequence Processing";		}
		
	string getHelpString();	
    string getOutputPattern(string);
    string getCommonQuestions();
		string getCitation() { return "http://www.mothur.org/wiki/Get.seqs"; }
		string getDescription()		{ return "gets sequences from a list, fasta, count, name, group, alignreport, quality, fastq, contigsreport or taxonomy file"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
        unordered_set <string> names;
        vector<string> fastafiles, namefiles, groupfiles, countfiles, alignfiles, listfiles, taxfiles, fastqfiles, contigsreportfiles, qualityfiles, outputNames;
		string accnosfile, accnosfile2, format, inputFileName;
		bool abort, dups;
        map<string, string> uniqueMap;
        map<string, set<string> > sanity; //for debug //maps file type to names chosen for file. something like "fasta" -> vector<string>. If running in debug mode this is filled and we check to make sure all the files have the same names. If they don't we output the differences for the user.
		
        void readFasta(map<string, vector<int> > nameToGroups, string fastafile, vector<string> outputFiles, vector<string>);
        void readFasta(string); //inputFastaFile, mothur generates output name
        void readFasta(string, string); //inputFastaFile, outputName (internal use)
        void readName(string); //inputNameFile, mothur generates output name
        void readName(string, string); //inputNameFile, outputName (internal use)
        void readCount(string); //inputCountFile, mothur generates output name
        void readCount(string, string); //inputCountFile, outputName (internal use)
        void readList(string); //inputListFile, mothur generates output name
        void readList(string, string); //inputListFile, outputName (internal use)
        void readList(string, vector<string>); //inputListFile, outputNames for each distance in list file (internal use)
        void readFastq(string);
        void readGZFastq(string);
        void readGroup(string);
        void readAlign(string);
        void readTax(string);
        void readQual(string);
        void readContigs(string);
		int compareAccnos(string);
        int runSanityCheck(string, string, string, string, string, string);
        int createMisMatchFile(ofstream&, string, string, set<string>, set<string>);
        int processList(ListVector*& list, string output, bool&);

		
};
//**********************************************************************************************************************

#endif


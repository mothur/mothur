#ifndef PARSEFASTAQCOMMAND_H
#define PARSEFASTAQCOMMAND_H

/*
 *  parsefastaqcommand.h
 *  Mothur
 *
 *  Created by westcott on 9/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "trimoligos.h"
#include "sequence.hpp"
#include "fastqread.h"
#include "groupmap.h"
#include "oligos.h"
#include "filefile.hpp"
#include "splitgroupscommand.h"


class ParseFastaQCommand : public Command {

public:
	ParseFastaQCommand(string);
	~ParseFastaQCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "fastq.info";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Parse.fastq"; }
	string getDescription()		{ return "reads a fastq file and creates a fasta and quality file"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:

	vector<string> outputNames;	
	string outputDir, inputDir, fastaQFile, format, oligosfile, groupfile, file, inputfile, ffastq, rfastq;
	bool abort, fasta, qual, pacbio, pairedOligos, reorient, createFileGroup, hasIndex;
    int pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, split, numBarcodes, numPrimers, numLinkers, numSpacers, numRPrimers, fileOption;
    GroupMap* groupMap;
    Oligos oligos;
    
    map<int, string> file2Group;
    vector< vector<string> > readFile();
    string ffqnoMatchFile, rfqnoMatchFile, ffnoMatchFile, rfnoMatchFile, fqnoMatchFile, rqnoMatchFile;
    vector<string> Groups;
    map<string, string> seqGroups;
    map<string, long long> groupCounts;

    set<string> processFile(string inputfile, TrimOligos*&, TrimOligos*&);
    int processFile(vector<string> inputfiles, TrimOligos*&, TrimOligos*&);
	vector<int> convertQual(string);
    vector<char> convertTable;
    bool readOligos(string oligosFile);
    bool readGroup(string oligosFile);
    int findGroup(Sequence&, QualityScores&, string&, TrimOligos*&, TrimOligos*&, int, int);
    int findGroup(Sequence, string&, string);
    int findGroup(Sequence&, QualityScores&, Sequence&, QualityScores&, string&, TrimOligos*&, TrimOligos*&, int, int);
    map<string, vector<string> > splitFastqFile(string outputGroupFile, string resultFastqfile); //uses split.groups command to parse the reads by sample
    
    
};

#endif



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
#include "groupmap.h"
#include "oligos.h"

struct fastqRead2 {
    string quality;
	Sequence seq;
    string wholeRead;
	
	fastqRead2() {  };
	fastqRead2(Sequence s, string q, string w) : seq(s), quality(q), wholeRead(w){};
	~fastqRead2() {};
};


class ParseFastaQCommand : public Command {

public:
	ParseFastaQCommand(string);
	ParseFastaQCommand();
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
	bool abort, fasta, qual, pacbio, pairedOligos, reorient, createFileGroup;
    int pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, split, numBarcodes, numPrimers, numLinkers, numSpacers, numRPrimers, fileOption;
    GroupMap* groupMap;
    Oligos oligos;
    
    map<int, string> file2Group;
    vector< vector<string> > readFile();
    vector<vector<string> > fastqFileNames;
    vector<vector<string> > rfastqFileNames;
    vector<vector<string> > fastaFileNames;
    vector<vector<string> > qualFileNames;
    vector<vector<string> > rfastaFileNames;
    vector<vector<string> > rqualFileNames;
    string ffqnoMatchFile, rfqnoMatchFile, ffnoMatchFile, rfnoMatchFile, fqnoMatchFile, rqnoMatchFile;
    vector<string> Groups;
    map<string, int> GroupToFile;

	
    int processFile(string inputfile, TrimOligos*&, TrimOligos*&);
    int processFile(vector<string> inputfiles, TrimOligos*&, TrimOligos*&);
	vector<int> convertQual(string);
    vector<char> convertTable;
    bool readOligos(string oligosFile);
    bool readGroup(string oligosFile);
    fastqRead2 readFastq(ifstream&, bool&);
    int findGroup(fastqRead2, int&, int&, TrimOligos*&, TrimOligos*&, int, int);
    int findGroup(fastqRead2, int&, int&, string);
    int findGroup(fastqRead2, fastqRead2, int&, int&, TrimOligos*&, TrimOligos*&, int, int);
    
    
};

#endif



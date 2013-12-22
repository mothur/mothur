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
	string outputDir, fastaQFile, format, oligosfile, groupfile;
	bool abort, fasta, qual, pacbio, pairedOligos;
    int pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, split;
    GroupMap* groupMap;
    
    //oligos file data structures
    vector<string> linker, spacer, primerNameVector, barcodeNameVector, revPrimer;
	map<string, int> barcodes;
	map<string, int> primers;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
    vector<vector<string> > fastqFileNames;
    string noMatchFile;
	
	vector<int> convertQual(string);
    vector<char> convertTable;
    bool readOligos(string oligosFile);
    bool readGroup(string oligosFile);
    string reverseOligo(string oligo);
    fastqRead2 readFastq(ifstream&, bool&);
    int findGroup(fastqRead2, int&, int&, TrimOligos*&, int, int);
    int findGroup(fastqRead2, int&, int&, string);
    
};

#endif



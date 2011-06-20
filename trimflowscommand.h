#ifndef TRIMFLOWSCOMMAND_H
#define TRIMFLOWSCOMMAND_H

/*
 *  trimflowscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"
#include "flowdata.h"
#include "groupmap.h"

class TrimFlowsCommand : public Command {
public:
	TrimFlowsCommand(string);
	TrimFlowsCommand();
	~TrimFlowsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "trim.flows";	}
	string getCommandCategory()		{ return "Hidden";		}
	string getHelpString();	
	string getCitation() { return "no citation"; }
	string getDescription()		{ return "trim.flows"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort;

	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};
	int comboStarts;
	vector<int> processIDS;   //processid
	vector<linePair*> lines;

	vector<unsigned long int> getFlowFileBreaks();
	int createProcessesCreateTrim(string, string, string, string, vector<vector<string> >); 
	int driverCreateTrim(string, string, string, string, vector<vector<string> >, linePair*);

	vector<string> outputNames;
	set<string> filesToRemove;
	
	void getOligos(vector<vector<string> >&);		//a rewrite of what is in trimseqscommand.h
	int stripBarcode(Sequence&, int&);				//largely redundant with trimseqscommand.h
	int stripForward(Sequence&, int&);				//largely redundant with trimseqscommand.h
	bool stripReverse(Sequence&);					//largely redundant with trimseqscommand.h
	bool compareDNASeq(string, string);				//largely redundant with trimseqscommand.h
	int countDiffs(string, string);					//largely redundant with trimseqscommand.h
	
	bool allFiles;
	int processors;
	int numFPrimers, numRPrimers;
	int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs;
	int numFlows;
	float signal, noise;
	bool fasta;
	string flowOrder;	
	
	string flowFileName, oligoFileName, outputDir;


	map<string, int> barcodes;
	map<string, int> primers;
	vector<string> revPrimer;

	vector<string> primerNameVector;	//needed here?
	vector<string> barcodeNameVector;	//needed here?

	map<string, int> combos;			//needed here?
	map<string, int> groupToIndex;		//needed here?
	
};


#endif

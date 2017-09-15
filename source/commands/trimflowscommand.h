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
#include "trimoligos.h"
#include "oligos.h"

class TrimFlowsCommand : public Command {
public:
	TrimFlowsCommand(string);
	TrimFlowsCommand();
	~TrimFlowsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "trim.flows";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Trim.flows"; }
	string getDescription()		{ return "trim.flows"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort;

	int comboStarts;
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
    vector<string> outputNames;
	set<string> filesToRemove;
    bool allFiles;
	int processors;
	int numFPrimers, numRPrimers, numBarcodes;
	int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs, sdiffs, ldiffs, numLinkers, numSpacers;
	int numFlows;
	float signal, noise;
	bool fasta, pairedOligos, reorient;
	string flowOrder, flowFileName, oligoFileName, outputDir;
    Oligos oligos;


	vector<unsigned long long> getFlowFileBreaks();
	int createProcessesCreateTrim(string, string, string, string, vector<vector<string> >); 
	int driverCreateTrim(string, string, string, string, vector<vector<string> >, linePair*);
	int getOligos(vector<vector<string> >&);		//a rewrite of what is in trimseqscommand.h
	
	
	
};
#endif

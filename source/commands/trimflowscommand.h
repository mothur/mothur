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
#include "splitgroupscommand.h"

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
    vector<string> outputNames;
	set<string> filesToRemove;
	int numFPrimers, numRPrimers, numBarcodes, processors, numFlows, comboStarts;
	int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs, sdiffs, ldiffs, numLinkers, numSpacers;
	float signal, noise;
	bool fasta, pairedOligos, reorient, allFiles, abort, createGroup;
	string flowOrder, flowFileName, oligoFileName, outputDir;
    
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
    map<string, int> barcodes;
    map<string, int> primers;
    vector<string>  linker;
    vector<string>  spacer;
    vector<string> primerNameVector;
    vector<string> barcodeNameVector;
    vector<string> revPrimer;
    map<string, string> groupMap;
    
	vector<double> getFlowFileBreaks();
	int createProcessesCreateTrim(string, string, string, string);
	int getOligos();
	
};
#endif

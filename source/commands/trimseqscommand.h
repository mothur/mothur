#ifndef TRIMSEQSCOMMAND_H
#define TRIMSEQSCOMMAND_H

/*
 *  trimseqscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"
#include "qualityscores.h"
#include "trimoligos.h"
#include "counttable.h"
#include "writer.h"
#include "splitgroupscommand.h"


class TrimSeqsCommand : public Command {
public:
	TrimSeqsCommand(string);
	TrimSeqsCommand();
	~TrimSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "trim.seqs";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Trim.seqs"; }
	string getDescription()		{ return "provides the preprocessing features needed to screen and sort pyrosequences"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:    
	bool abort, createGroup;
	string fastaFile, oligoFile, qFileName, groupfile, nameFile, countfile, outputDir;
	
	bool flip, allFiles, qtrim, keepforward, pairedOligos, reorient, logtransform;
	int numFPrimers, numRPrimers, numLinkers, numSpacers, maxAmbig, maxHomoP, minLength, maxLength, processors, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs, comboStarts;
	int qWindowSize, qWindowStep, keepFirst, removeLast;
	double qRollAverage, qThreshold, qWindowAverage, qAverage;
	vector<string> revPrimer, outputNames;
	set<string> filesToRemove;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
	map<string, int> barcodes;
	vector<string> groupVector;
	map<string, int> primers;
    vector<string>  linker;
    vector<string>  spacer;
	map<string, int> combos;
	map<string, int> groupToIndex;
	vector<string> primerNameVector;	//needed here?
	vector<string> barcodeNameVector;	//needed here?
	map<string, int> groupCounts;  
	map<string, string> nameMap;
    map<string, int> nameCount; //for countfile name -> repCount
    map<string, string> groupMap; //for countfile name -> group

	vector<int> processIDS;   //processid
	vector<linePair> lines;
	vector<linePair> qLines;
	
    bool getOligos();
	int createProcessesCreateTrim(string, string, string, string, string, string, string, string, string, string, string);
    int processNamesCountFiles(string trimFasta, set<string> badNames, map<string, string> groupMap, string trimNameFileName, string scrapNameFileName, string trimCountFileName, string scrapCountFileName, string groupFile);
	int setLines(string, string);
};

/**************************************************************************************************/

#endif


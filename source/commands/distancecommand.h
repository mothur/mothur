#ifndef DISTANCECOMMAND_H
#define DISTANCECOMMAND_H

/*
 *  distancecommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "validcalculator.h"
#include "calculator.h"
#include "sequencedb.h"
#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"
#include "jtt.hpp"
#include "writer.h"
#include "proteindb.hpp"
#include "pmb.hpp"
#include "pam.hpp"
#include "kimura.hpp"

/**************************************************************************************************/
struct distanceData {
    long long startLine, endLine, numNewFasta, count;
	float cutoff;
    StorageDatabase* db;
    StorageDatabase* oldFastaDB;
	MothurOut* m;
	OutputWriter* threadWriter;
    string outputFileName, calc;
	bool countends, prot;
    Utils util;
	
	distanceData(){}
    distanceData(OutputWriter* ofn) {
        threadWriter = ofn;
        m = MothurOut::getInstance();
    }
    
    distanceData(string ofn) {
        outputFileName = ofn;
        m = MothurOut::getInstance();
    }
	void setVariables(int s, int e,  float c, StorageDatabase*& dbsp, StorageDatabase*& oldfn, string Est, bool met, long long num, bool cnt) {
		startLine = s;
		endLine = e;
		cutoff = c;
		db = dbsp;
        oldFastaDB = oldfn;
		calc = Est;
        prot = met;
		numNewFasta = num;
		countends = cnt;
        count = 0;
	}
};
/**************************************************************************************************/
class DistanceCommand : public Command {

public:
    DistanceCommand(string);
    DistanceCommand(StorageDatabase*&, string, double, string, int); //used by mothur's splitMatrix class to avoid rereading files
	~DistanceCommand() = default;
	
	vector<string> setParameters();
	string getCommandName()			{ return "dist.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD (2010). The effects of alignment quality, distance calculation method, sequence filtering, and region on the analysis of 16S rRNA gene-based studies. PLoS Comput Biol 6: e1000844. \nhttp://www.mothur.org/wiki/Dist.seqs"; }
	string getDescription()		{ return "calculate the pairwaise distances between aligned sequences"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	
    StorageDatabase* db;
	string output, fastafile, calc,  oldfastafile, column;
    int processors;
    long long numNewFasta, numSeqs, numDistsBelowCutoff;
	float cutoff;
	
	bool abort, countends, fitCalc, prot, compress;
	vector<string> outputNames; 
	
	void createProcesses(string);
	bool sanityCheck();
};

#endif

/**************************************************************************************************/




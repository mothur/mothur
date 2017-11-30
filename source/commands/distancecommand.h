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

//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct distanceData {
	long long startLine, endLine;
	string dFileName;
	float cutoff;
	SequenceDB alignDB;
	vector<string> Estimators;
	MothurOut* m;
	string output;
	long long numNewFasta, count;
	bool countends;
    Utils util;
	
	distanceData(){}
	distanceData(int s, int e, string dbname, float c, SequenceDB db, vector<string> Est, MothurOut* mout, string o, long long num, bool cnt) {
		startLine = s;
		endLine = e;
		dFileName = dbname;
		cutoff = c;
		alignDB = db;
		Estimators = Est;
		m = mout;
		output = o;
		numNewFasta = num;
		countends = cnt;
        count = 0;
		
	}
};
/**************************************************************************************************/
class DistanceCommand : public Command {

public:
	DistanceCommand(string);
	DistanceCommand();
	~DistanceCommand() {}
	
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
	
    SequenceDB alignDB;
	string output, fastafile, calc, outputDir, oldfastafile, column, compress;
    int processors;
    long long numNewFasta, numSeqs;
	float cutoff;
	
	bool abort, countends;
	vector<string>  Estimators, outputNames; //holds estimators to be used
	
	void createProcesses(string);
	bool sanityCheck();
};

#endif

/**************************************************************************************************/




#ifndef SENSSPECCOMMAND_H
#define SENSSPECCOMMAND_H


/*
 *  sensspeccommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 7/6/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "listvector.hpp"
#include "inputdata.h"

class SensSpecCommand : public Command {
	
public:
	SensSpecCommand(string);
	SensSpecCommand();
	~SensSpecCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "sens.spec";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "Schloss PD, Westcott SL (2011). Assessing and improving methods used in OTU-based approaches for 16S rRNA gene sequence analysis. Appl Environ Microbiol. \nhttp://www.mothur.org/wiki/Sens.spec"; }
	string getDescription()		{ return "sens.spec"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	int processPhylip();
	int processColumn();
	void setUpOutput();
	void outputStatistics(string, string);
	
	string listFile, distFile, nameFile, sensSpecFileName, phylipfile, columnfile;
	string outputDir;
	string format;
	vector<string> outputNames;
	set<string> labels; //holds labels to be used

	long int truePositives, falsePositives, trueNegatives, falseNegatives;
	bool abort, allLines;
	bool hard;
	//string lineLabel;
	double cutoff;
	int precision;
	
	int fillSeqMap(map<string, int>&, ListVector*&);
	int fillSeqPairSet(set<string>&, ListVector*&);
	int process(map<string, int>&, string, bool&, string&);
	int process(set<string>&, string, bool&, string&, int);

};

#endif




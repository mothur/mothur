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

class SensSpecCommand : public Command {
	
public:
	SensSpecCommand(string);
	SensSpecCommand();
	~SensSpecCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "sens.spec";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Sens.spec"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	void processPhylip();
	void processColumn();
	void setUpOutput();
	void outputStatistics(string, string);
	
	string listFile, distFile, nameFile, sensSpecFileName, phylipfile, columnfile;
	string outputDir;
	string format;
	vector<string> outputNames;

	long int truePositives, falsePositives, trueNegatives, falseNegatives;
	bool abort;
	bool hard;
	string lineLabel;
	double cutoff;
	int precision;
};

#endif




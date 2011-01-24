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
	~SensSpecCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	void processPhylip();
	void processColumn();
	void setUpOutput();
	void outputStatistics(string, string);
	
	string listFile, distFile, nameFile, sensSpecFileName;
	string outputDir;
	string format;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;

	long int truePositives, falsePositives, trueNegatives, falseNegatives;
	bool abort;
	bool hard;
	string lineLabel;
	double cutoff;
	int precision;
};

#endif




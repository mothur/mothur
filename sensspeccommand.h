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
	~SensSpecCommand();
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

	long int truePositives, falsePositives, trueNegatives, falseNegatives;
	bool abort;
	bool hard;
	string lineLabel;
	double cutoff;
	int precision;
};

#endif
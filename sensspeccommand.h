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
///	void processColumn(map<string, int>);
	
	string listFile, distFile, nameFile, outputDir;
	string format;
//	int numSeqs, numDists;
	int truePositives, falsePositives, trueNegatives, falseNegatives;
	bool abort;
	bool hard;
	string lineLabel;
	double cutoff;
};

#endif
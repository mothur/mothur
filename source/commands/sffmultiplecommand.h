#ifndef Mothur_sffmultiplecommand_h
#define Mothur_sffmultiplecommand_h

//
//  sffmultiplecommand.h
//  Mothur
//
//  Created by Sarah Westcott on 8/14/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "command.hpp"
#include "sffinfocommand.h"
#include "seqsummarycommand.h"
#include "trimflowscommand.h"
#include "shhhercommand.h"
#include "trimseqscommand.h"

class SffMultipleCommand : public Command {
	
public:
	SffMultipleCommand(string);
	~SffMultipleCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "sff.multiple";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Sff.multiple"; }
	string getDescription()		{ return "run multiple sff files through, sffinfo, trim.flow, shhh.flows and trim.seqs combining the results"; }
    
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    string inputDir;
	string filename,  flowOrder, lookupFileName, minDelta;
	vector<string> outputNames;
	bool abort, trim, large, flip, allFiles, keepforward, append, makeGroup;
	int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs, sdiffs, ldiffs;
	int processors, maxIters, largeSize;
	float signal, noise, cutoff, sigma;
    int keepFirst, removeLast, maxAmbig;
    
    int readFile(vector<string>& sffFiles, vector<string>& oligosFiles);
    long long createProcesses(vector<string> sffFiles, vector<string> oligosFiles, string, string, string);
};

#endif

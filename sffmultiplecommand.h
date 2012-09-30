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
	SffMultipleCommand();
	~SffMultipleCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "sff.multiple";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Sff.multiple"; }
	string getDescription()		{ return "run multiple sff files through, sffinfo, trim.flow, shhh.flows and trim.seqs combining the results"; }
    
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    
    struct linePair {
		int start;
		int end;
		linePair(int i, int j) : start(i), end(j) {}
	};

	string filename, outputDir, flowOrder, lookupFileName, minDelta;
	vector<string> outputNames;
	bool abort, trim, large, flip, allFiles, keepforward, append, makeGroup;
	int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs, sdiffs, ldiffs;
	int processors, maxIters, largeSize;
	float signal, noise, cutoff, sigma;
    int keepFirst, removeLast, maxAmbig;
    
    int readFile(vector<string>& sffFiles, vector<string>& oligosFiles);
    int createProcesses(vector<string> sffFiles, vector<string> oligosFiles, string, string, string);
    int driver(vector<string> sffFiles, vector<string> oligosFiles, int start, int end, string, string, string);
    int mergeOutputFileList(map<string, vector<string> >& files, map<string, vector<string> >& temp);


    
};

#endif

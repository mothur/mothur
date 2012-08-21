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

	string filename, outputDir, flowOrder;
	vector<string> outputNames;
	bool abort, trim, large, flip, qtrim, allFiles, keepforward;
	int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs, sdiffs, ldiffs, numLinkers, numSpacers;
	int numFlows, numFPrimers, numRPrimers, processors, maxIters, largeSize;
	float signal, noise, cutoff, sigma, minDelta;
    int qWindowSize, qWindowStep, keepFirst, removeLast, maxAmbig;
	double qRollAverage, qThreshold, qWindowAverage, qAverage;
    
    int readFile(vector<string>& sffFiles, vector<string>& oligosFiles);
    int createProcesses(vector<string> sffFiles, vector<string> oligosFiles);
    int driver(vector<string> sffFiles, vector<string> oligosFiles, int start, int end);

    
};


#endif

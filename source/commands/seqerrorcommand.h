#ifndef SEQERRORCOMMAND
#define SEQERRORCOMMAND

/*
 *  seqerrorcommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 7/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "sequence.hpp"
#include "counttable.h"
#include "compare.h"


class SeqErrorCommand : public Command {
public:
	SeqErrorCommand(string);
	SeqErrorCommand();
	~SeqErrorCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "seq.error";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Gevers D, Westcott SL (2011).  Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies.  PLoS ONE.  6:e27310.\nhttp://www.mothur.org/wiki/Seq.error"; }
	string getDescription()		{ return "seq.error"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort;
	
	vector<Sequence> getReferences(string);
	void printSubMatrix(vector<vector<int> >& substitutionMatrix);
	void printErrorFRFile(map<char, vector<int> >& errorForward, map<char, vector<int> >& errorReverse);
	void printErrorQuality(map<char, vector<int> >&);
	void printQualityFR(vector<vector<int> >& qualForwardMap, vector<vector<int> >& qualReverseMap);
	
	int setLines(string, string, string, vector<linePair>&, vector<linePair>&, vector<linePair>&);
    int createProcesses(string, string, string, string, string, string, vector<vector<int> >&, vector<vector<int> >&, vector<vector<int> >&, vector<int>&,  map<char, vector<int> >&, map<char, vector<int> >&, map<char, vector<int> >&, vector<string>&, vector<Sequence>&);

	string queryFileName, referenceFileName, qualFileName, reportFileName, namesFileName, outputDir, countfile;
	double threshold;
	bool ignoreChimeras, aligned;
	int numRefs, processors;
	int maxLength, totalBases, totalMatches;
	vector<string> outputNames;
};

#endif

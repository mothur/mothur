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
	
	vector<int> processIDS;   //processid
	vector<linePair> lines;
	vector<linePair> qLines;
	vector<linePair> rLines;

	void getReferences();
	map<string,int> getWeights();
	Compare getErrors(Sequence, Sequence);
	void printErrorHeader(ofstream&);
	void printErrorData(Compare, int, ofstream&, ofstream&);
	void printSubMatrix();
	void printErrorFRFile(map<char, vector<int> >, map<char, vector<int> >);
	void printErrorQuality(map<char, vector<int> >);
	void printQualityFR(vector<vector<int> >, vector<vector<int> >);
	
	int setLines(string, string, string);
	int driver(string, string, string, string, string, string, linePair, linePair, linePair);
	int createProcesses(string, string, string, string, string, string);

	string queryFileName, referenceFileName, qualFileName, reportFileName, namesFileName, outputDir, countfile;
	double threshold;
	bool ignoreChimeras, aligned;
	int numRefs, processors;
	int maxLength, totalBases, totalMatches;
	//ofstream errorSummaryFile, errorSeqFile;
	vector<string> outputNames;
	
	vector<Sequence> referenceSeqs;
	vector<vector<int> > substitutionMatrix;
	vector<vector<int> > qualForwardMap;
	vector<vector<int> > qualReverseMap;
	vector<int> misMatchCounts;
	map<char, vector<int> > qScoreErrorMap;
	map<char, vector<int> > errorForward;
	map<char, vector<int> > errorReverse;
	map<string, int> weights;
	vector<string> megaAlignVector;

};

#endif

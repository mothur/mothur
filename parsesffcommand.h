#ifndef PARSESFFCOMMAND_H
#define PARSESFFCOMMAND_H

/*
 *  parsesffcommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 2/6/10.
 *  Copyright 2010 Patrick D. Schloss. All rights reserved.
 *
 */


#include "command.hpp"

class ParseSFFCommand : public Command {
public:
	ParseSFFCommand(string);
	~ParseSFFCommand();
	int execute();
	void help();	
	
private:

	int parseHeaderLineToInt(ifstream&);
	vector<float> parseHeaderLineToFloatVector(ifstream&, int);
	vector<int> parseHeaderLineToIntVector(ifstream&, int);
	string parseHeaderLineToString(ifstream&);
	void screenFlow(vector<float>, int&);
	string flow2seq(vector<float>, int);
	bool screenSeq(string&, int&);
	bool compareDNASeq(string, string);
	void getOligos(vector<ofstream*>&);
	
	
	string sffFile;
	string oligoFile;

	int minLength;
	int numFPrimers, numRPrimers, numBarcodes;
	vector<string> forPrimer, revPrimer;
	map<string, int> barcodes;
	vector<string> groupVector;
	vector<string> outputNames;

//	string stripSeqQual(string, int, int);
//	string stripQualQual(string, int, int);
	
	string outputDir;
	bool abort;
};

#endif



#ifndef ANOSIMCOMMAND_H
#define ANOSIMCOMMAND_H

/*
 *  anosimcommand.h
 *  mothur
 *
 *  Created by westcott on 2/14/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "inputdata.h"
#include "sharedrabundvector.h"
#include "validcalculator.h"
#include "readphylipvector.h"

class GlobalData;

class AnosimCommand : public Command {
	
public:
	AnosimCommand(string);
	AnosimCommand();
	~AnosimCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	struct linePair {
		int start;
		int num;
		linePair(int i, int j) : start(i), num(j) {}
	};
	vector<linePair> lines;
	
	GlobalData* globaldata;
	GroupMap* designMap;
	map<string, vector<string> > outputTypes;
	
	vector< vector<double> > matrix;
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string format, groups, label, outputDir, inputDir, designfile, sets, phylipfile, calc, sharedfile;
	vector<string> Groups, outputNames, Sets;
	vector< vector<string> > namesOfGroupCombos;
	int iters, processors;
	vector<Calculator*> calculators;
	
	int driver(int, int, vector<SharedRAbundVector*>, string);
	int driver(int, int, vector<string>, string, vector< vector<double> >&);
	int process(vector<SharedRAbundVector*>);
	int calcAnosim(ofstream&, int, vector<string>);
	double calcWithin(vector<seqDist>&, vector<string>);
	double calcBetween(vector<seqDist>&);
	vector<seqDist> convertToRanks();
};

#endif




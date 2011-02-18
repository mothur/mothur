#ifndef AMOVACOMMAND_H
#define AMOVACOMMAND_H

/*
 *  amovacommand.h
 *  mothur
 *
 *  Created by westcott on 2/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "sharedrabundvector.h"
#include "validcalculator.h"
#include "readphylipvector.h"

class GlobalData;

class AmovaCommand : public Command {
	
public:
	AmovaCommand(string);
	AmovaCommand();
	~AmovaCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	double runAMOVA(ofstream&, map<string, vector<int> >, double);
	double calcSSWithin(map<string, vector<int> >&);
	double calcSSTotal(map<string, vector<int> >&);
	map<string, vector<int> > getRandomizedGroups(map<string, vector<int> >);

	
	bool abort;
	GlobalData* globaldata;
	map<string, vector<string> > outputTypes;
	vector<string> outputNames;

	string outputDir, inputDir, designFileName, phylipFileName;
	GroupMap* designMap;
	vector< vector<double> > distanceMatrix;
	int iters;
	double experimentwiseAlpha;
	
//	struct linePair {
//		int start;
//		int num;
//		linePair(int i, int j) : start(i), num(j) {}
//	};
//	vector<linePair> lines;
//
//	vector< vector<string> > namesOfGroupCombos;
//	vector<string> Groups, outputNames, Sets;
//	int processors;
//	string groups, sets, calc, sharedfile, label, allLines, pickedGroups;
//	vector<Calculator*> calculators;
//	set<string> labels; //holds labels to be used
//	int driver(int, int, vector<SharedRAbundVector*>, string);
//	int driver(int, int, vector<string>, string, vector< vector<double> >&);
//	int process(vector<SharedRAbundVector*>);	
};

#endif


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

class GroupMap;

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
	bool abort;
	GroupMap* designMap;
	map<string, vector<string> > outputTypes;
	string outputDir, inputDir, designFileName, phylipFileName;
	
	vector<vector<double> > convertToRanks(vector<vector<double> >);
	double calcR(vector<vector<double> >, map<string, vector<int> >);
	map<string, vector<int> > getRandomizedGroups(map<string, vector<int> >);
	double runANOSIM(ofstream&, vector<vector<double> >, map<string, vector<int> >, double);
	
	vector< vector<double> > distanceMatrix;
	vector<string> outputNames;
	int iters;
	double experimentwiseAlpha;
	vector< vector<string> > namesOfGroupCombos;
	
	
};

#endif




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

//class GlobalData;
class GroupMap;

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
	map<string, vector<string> > outputTypes;
	vector<string> outputNames;

	string outputDir, inputDir, designFileName, phylipFileName;
	GroupMap* designMap;
	vector< vector<double> > distanceMatrix;
	int iters;
	double experimentwiseAlpha;
};

#endif


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
	~AnosimCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "anosim";					}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getHelpString();	
	string getCitation() { return "Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in community structure. _Australian Journal of Ecology_ 18, 117-143.\nhttp://www.mothur.org/wiki/Anosim"; }
	string getDescription()		{ return "anosim"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:
	bool abort;
	GroupMap* designMap;
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




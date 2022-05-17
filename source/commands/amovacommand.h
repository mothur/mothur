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
class DesignMap;

class AmovaCommand : public Command {
	
public:
	AmovaCommand(string);
	~AmovaCommand() = default;
	
	vector<string> setParameters();
	string getCommandName()			{ return "amova";					}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();
    string getCommonQuestions();
    string getOutputPattern(string);	
    string getCitation() { return "Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecol 26: 32-46.\nhttp://www.mothur.org/wiki/Amova"; }
	string getDescription()		{ return "analysis of molecular variance"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:
	double runAMOVA(ofstream&, map<string, vector<int> >, double);
	double calcSSWithin(map<string, vector<int> >&);
	double calcSSTotal(map<string, vector<int> >&);
	map<string, vector<int> > getRandomizedGroups(map<string, vector<int> >);

	bool abort;
	vector<string> outputNames, Sets;

	string inputDir, designFileName, phylipFileName;
	DesignMap* designMap;
	vector< vector<double> > distanceMatrix;
	int iters;
	double experimentwiseAlpha;
};

#endif


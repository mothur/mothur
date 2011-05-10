#ifndef MANTELCOMMAND_H
#define MANTELCOMMAND_H

/*
 *  mantelcommand.h
 *  mothur
 *
 *  Created by westcott on 2/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "linearalgebra.h"

class MantelCommand : public Command {
public:
	MantelCommand(string);
	MantelCommand();
	~MantelCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "mantel";					}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getHelpString();	
	string getCitation() { return "Sokal, R. R., & Rohlf, F. J. (1995). Biometry, 3rd edn. New York: Freeman. http://www.mothur.org/wiki/Mantel"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:
	
	string phylipfile1, phylipfile2, outputDir, method;
	bool abort;
	int iters;
	
	vector<string> outputNames;
};


#endif




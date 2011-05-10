#ifndef PCOACOMMAND_H
#define PCOACOMMAND_H

/*
 *  pcoacommand.h
 *  Mothur
 *
 *  Created by westcott on 1/4/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "linearalgebra.h"


/*****************************************************************/
class PCOACommand : public Command {
	
public:
	PCOACommand(string);	
	PCOACommand();
	~PCOACommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pcoa";					}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Pcoa"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:

	bool abort, metric;
	string phylipfile, filename, fbase, outputDir;
	vector<string> outputNames;
	LinearAlgebra linearCalc;
	
	void get_comment(istream&, char, char);
	void output(string, vector<string>, vector<vector<double> >&, vector<double>);
	
};
	
/*****************************************************************/
	
#endif


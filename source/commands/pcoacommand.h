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
	~PCOACommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pcoa";					}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "McCune B, Grace JB, Urban DL (2002). Analysis of ecological communities. MjM Software Design: Gleneden Beach, OR. \nLegendre P, Legendre L (1998). Numerical Ecology. Elsevier: New York. \nhttp://www.mothur.org/wiki/Pcoa"; }
	string getDescription()		{ return "pcoa"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:

	bool abort, metric;
	string phylipfile, filename, fbase;
	vector<string> outputNames;
	LinearAlgebra linearCalc;
	
	void get_comment(istream&, char, char);
	void output(string, vector<string>, vector<vector<double> >&, vector<double>);
	
};
	
/*****************************************************************/
	
#endif


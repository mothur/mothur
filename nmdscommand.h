#ifndef NMDSCOMMAND_H
#define NMDSCOMMAND_H

/*
 *  nmdscommand.h
 *  mothur
 *
 *  Created by westcott on 1/11/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "linearalgebra.h"


/*  
 Translated from the nmds.R code written by Sarah Goslee using,
 
 # Non-metric multidimensional scaling function
 # using the majorization algorithm from
 # Borg & Groenen 1997, Modern Multidimensional Scaling.
 #
 
 # also referenced (Kruskal 1964)
 
 */

/*****************************************************************/
class NMDSCommand : public Command {
	
public:
	NMDSCommand(string);	
	NMDSCommand();
	~NMDSCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "nmds";					}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getHelpString();	
	string getCitation() { return "Borg, Groenen (1997). Non-metric multidimensional scaling function using the majorization algorithm, in Modern Multidimensional Scaling. Ed. T.F. Cox and M.A.A. Cox. Chapman and Hall. \nhttp://www.mothur.org/wiki/Nmds"; }
	string getDescription()		{ return "nmds"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:
	
	bool abort;
	string phylipfile, outputDir, axesfile;
	int maxdim, mindim, maxIters, iters;
	double epsilon;
	vector<string> outputNames;
	LinearAlgebra linearCalc;
	
	vector< vector<double> > nmdsCalc(vector< vector<double> >&, vector< vector<double> >&, double&);
	vector< vector<double> > getConfiguration(vector< vector<double> >&, int);
	vector< vector<double> > generateStartingConfiguration(int, int); //pass in numNames, return axes
	int normalizeConfiguration(vector< vector<double> >&, int, int);
	double calculateStress(vector< vector<double> >&, vector< vector<double> >&);
	vector< vector<double> > readAxes(vector<string>);
	int output(vector< vector<double> >&, vector<string>&, ofstream&);	
};

/*****************************************************************/

#endif



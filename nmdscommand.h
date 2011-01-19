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
	~NMDSCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();	
	void help();
	
private:
	
	bool abort, trace;
	string phylipfile, outputDir, axesfile;
	int maxdim, mindim, maxIters, iters;
	double epsilon;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	LinearAlgebra linearCalc;
	
	vector< vector<double> > nmdsCalc(vector< vector<double> >&, vector< vector<double> >&, double&);
	vector< vector<double> > getConfiguration(vector< vector<double> >&, int);
	vector< vector<double> > generateStartingConfiguration(int, int); //pass in numNames, return axes
	int normalizeConfiguration(vector< vector<double> >&, int, int);
	double calculateStress(vector< vector<double> >&, vector< vector<double> >&);
	vector< vector<double> > readAxes(vector<string>);
	int output(vector< vector<double> >&, vector<string>&, ofstream&);

	//vector<seqDist> satisfyMonotonicity(vector<seqDist>, vector<int>);
	//vector< vector<double> > calculateStressGradientVector(vector<seqDist>&, vector<seqDist>&, double, double, vector< vector<double> >&);
	//double calculateMagnitude(vector< vector<double> >&);
	//double calculateStep(vector< vector<double> >&, vector< vector<double> >&, vector<double>&);
	//vector< vector<double> > calculateNewConfiguration(double, vector< vector<double> >&, vector< vector<double> >&);
	
};

/*****************************************************************/

#endif



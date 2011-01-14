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


/*  references used to make this command: "Nonmetric Multidimensional Scalling: A Numerical Method"
 by J. B. Kruskal  Psychometrika - Vol 29, No. 2 June 1964 */

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
	
	bool abort;
	string phylipfile, outputDir, axesfile;
	int dimension, maxIters;
	double step, cutoff;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	LinearAlgebra linearCalc;
	
	vector< vector<double> > generateStartingConfiguration(int); //pass in numNames, return axes
	int normalizeConfiguration(vector< vector<double> >&, int);
	vector<seqDist> satisfyMonotonicity(vector<seqDist>);
	double calculateStress(vector<seqDist>&, vector<seqDist>&, double&);
	vector< vector<double> > calculateStressGradientVector(vector<seqDist>&, vector<seqDist>&, double, double, vector< vector<double> >&);
	double calculateMagnitude(vector< vector<double> >&);
	double calculateStep(vector< vector<double> >&, vector< vector<double> >&, vector<double>&);
	vector< vector<double> > calculateNewConfiguration(double, vector< vector<double> >&, vector< vector<double> >&);
	vector< vector<double> > readAxes(vector<string>);
	int output(string, string, vector< vector<double> >&, vector<double>&, vector<string>&);
};

/*****************************************************************/

#endif



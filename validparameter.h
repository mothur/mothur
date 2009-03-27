#ifndef VALIDPARAMETERS_H
#define VALIDPARAMETERS_H

/*
 *  validparameter.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
using namespace std;

#include "mothur.h"
#include "utilities.hpp"

//This class contains a list of all valid parameters in Mothur.  
//It has a function which will tell you if your parameter is valid.
//When adding a new parameter you must add it to the valid list in the class constructor.


class ValidParameters {

	public:
		ValidParameters();
		~ValidParameters();
		bool isValidParameter(string);
		bool isValidParameter(string, string, string);
		vector <string> addParameters(string[], int);
		void initCommandParameters();
		void initParameterRanges();

	private:
		map<string, string> readdist;
		map<string, string> readotu;
		map<string, string> readtree;
		map<string, string> cluster;
		map<string, string> deconvolute;
		map<string, string> parsimony;
		map<string, string> collectsingle;
		map<string, string> collectshared;
		map<string, string> rarefactsingle;
		map<string, string> rarefactshared;
		map<string, string> summarysingle;
		map<string, string> summaryshared;
		map<string, string> unifracweighted;
		map<string, string> unifracunweighted;
		map<string, string> libshuff;
		map<string, string> heatmap;
		
		map<string, string>::iterator it;
		map<string, vector<string> > commandParameters;
		map<string, vector<string> > parameterRanges;

};

#endif

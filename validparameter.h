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

//This class contains a list of all valid parameters in Mothur.  
//It has a function which will tell you if your parameter is valid.
//When adding a new parameter you must add it to the valid list in the class constructor.


class ValidParameters {

	public:
		ValidParameters();
		~ValidParameters();
		bool isValidParameter(string, string);
		
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
		
		void initialReaddist();
		void initialReadotu();
		void initialReadtree();
		void initialCluster();
		void initialDeconvolute();
		void initialParsimony();
		void initialCollectsingle();
		void initialCollectshared();
		void initialRarefactsingle();
		void initialRarefactshared();
		void initialSummarysingle();
		void initialSummaryshared();
		void initialUnifracweighted();
		void initialUnifracunweighted();
		void initialLibshuff();
		void initialHeatmap();

};

#endif

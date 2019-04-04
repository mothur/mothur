#ifndef VALIDCALCULATOR_H
#define VALIDCALCULATOR_H

/*
 *  validcalculator.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"

//This class contains a list of all valid calculators in Mothur.  
//It has a function which will tell you if your calculator is valid for the given parameter.
//When adding a new calculator you must add it to the valid list.


class ValidCalculators {
	public:
		ValidCalculators();
		~ValidCalculators();
		bool isValidCalculator(string, string);
		void printCalc(string, ostream&);
		string printCalc(string);
		void printCitations(vector<string>);
		
	private:
        map<string, string> estimators;
		map<string, string> single;
		map<string, string> shared;
		map<string, string> rarefaction;
		map<string, string> summary;
		map<string, string> sharedrarefaction;
		map<string, string> sharedsummary;
		map<string, string> vennsingle;
		map<string, string> vennshared;
		map<string, string> treegroup;
		map<string, string> matrix;
		map<string, string> heat;
		map<string, string> boot;
		map<string, string> distance;
		map<string, string>::iterator it;
		set<string> allCalcs;
		
		void initialSingle();
		void initialShared();
		void initialRarefaction();
		void initialSharedRarefact();
		void initialSummary();
		void initialSharedSummary();
		void initialVennSingle();
		void initialVennShared();
		void initialTreeGroups();
		void initialMatrix();
		void initialBoot();
		void initialDistance();
		void initialHeat();
        void initialEstimators();
    
		MothurOut* m;
};

#endif

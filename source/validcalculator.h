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

#include "mothurout.h"
#include "utils.hpp"

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
        set<string> estimators;
		set<string> single;
		set<string> shared;
		set<string> rarefaction;
		set<string> summary;
		set<string> sharedrarefaction;
		set<string> sharedsummary;
		set<string> vennsingle;
		set<string> vennshared;
		set<string> treegroup;
		set<string> matrix;
        set<string> clr;
		set<string> heat;
		set<string> distance;
        set<string> protdistance;
		set<string>::iterator it;
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
        void initialCLR();
		void initialDistance();
        void initialProtDistance();
		void initialHeat();
        void initialEstimators();
    
		MothurOut* m;
};

#endif

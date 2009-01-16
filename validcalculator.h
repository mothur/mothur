/*
 *  validcalculator.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
using namespace std;

#include <Carbon/Carbon.h>
#include <string>
#include <iostream>
#include <map>

//This class contains a list of all valid calculators in Mothur.  
//It has a function which will tell you if your calculator is valid for the given parameter.
//When adding a new calculator you must add it to the valid list.


class ValidCalculators {
	public:
		ValidCalculators();
		~ValidCalculators();
		bool isValidCalculator(string, string);
		
	private:
		map<string, string> single;
		map<string, string> shared;
		map<string, string> rarefaction;
		map<string, string> summary;
		map<string, string> sharedrarefaction;
		map<string, string> sharedsummary;
		void initialSingle();
		void initialShared();
		void initialRarefaction();
		void initialSharedRarefact();
		void initialSummary();
		void initialSharedSummary();
		
		
};

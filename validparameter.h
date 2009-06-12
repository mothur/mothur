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

#include "mothur.h"

//This class contains a list of all valid parameters in Mothur.  
//It has a function which will tell you if your parameter is valid.
//When adding a new parameter you must add it to the valid list in the class constructor.


class ValidParameters {

	public:
		ValidParameters();
		~ValidParameters();
		//bool isValidParameter(string, string, string) {return true;}
		bool isValidParameter(string, vector<string>, string);
		vector <string> addParameters(string[], int);
		void initParameterRanges();
		string validFile(map<string, string>, string, bool); //container, parameter, isFile

	private:
		map<string, string>::iterator it;
		map<string, vector<string> > parameterRanges;

};

#endif

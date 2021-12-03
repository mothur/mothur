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

#include "currentfile.h"
#include "mothurout.h"

//This class contains a list of all valid parameters in Mothur.  
//It has a function which will tell you if your parameter is valid.
//When adding a new parameter you must add it to the valid list in the class constructor.


class ValidParameters {

	public:
		ValidParameters();
		ValidParameters(string);
		~ValidParameters();
		bool isValidParameter(string, vector<string>, string);
		vector <string> addParameters(string[], int);
		void initParameterRanges();
		vector<string> validFiles(map<string, string>&, string);
        vector<string> validFastqGZFiles(map<string, string>&, string, bool&); 
        string validFile(map<string, string>&, string);
        string valid(map<string, string>&, string);
        string validPath(map<string, string>&, string);

	private:
		map<string, string>::iterator it;
		map<string, vector<string> > parameterRanges;
		MothurOut* m;
        CurrentFile* current;
        vector< vector<string> > locations;

};

#endif

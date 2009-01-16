/*
 *  validparameter.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
using namespace std;

#include <Carbon/Carbon.h>
#include <iostream>
#include <string>
#include <map>

//This class contains a list of all valid parameters in Mothur.  
//It has a function which will tell you if your parameter is valid.
//When adding a new parameter you must add it to the valid list in the class constructor.


class ValidParameters {

	public:
		ValidParameters();
		~ValidParameters();
		bool isValidParameter(string);
	private:
		map<string, string> parameters;

};
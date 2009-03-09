/*
 *  validparameter.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "validparameter.h"

/***********************************************************************/

ValidParameters::ValidParameters() {
	try {
	
		parameters["phylip"]		= "phylip";
		parameters["column"]		= "column";
		parameters["list"]			= "list"; 
		parameters["rabund"]		= "rabund"; 
		parameters["sabund"]		= "sabund";
		parameters["shared"]		= "shared";
		parameters["name"]			= "name"; 
		parameters["group"]			= "group"; 
		parameters["order"]			= "order"; 
		parameters["fasta"]			= "fasta"; 
		parameters["tree"]			= "tree";
		parameters["fileroot"]			= "fileroot";
		parameters["cutoff"]			= "cutoff"; 
		parameters["method"]			= "method";
		parameters["format"]			= "format"; 
		parameters["precision"]			= "precision"; 
		parameters["label"]				= "label"; 
		parameters["line"]				= "line";
		parameters["iters"]				= "iters"; 
		parameters["jumble"]			= "jumble"; 
		parameters["freq"]				= "freq"; 
		parameters["abund"]             = "abund";
		parameters["random"]			= "random";
		parameters["groups"]			= "groups";
		parameters["calc"]				= "calc";

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function ValidParameters. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function ValidParameters. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

ValidParameters::~ValidParameters() {}

/***********************************************************************/
bool ValidParameters::isValidParameter(string parameter) {
	try {	
	
		//is the parameter in the map
		if ((parameters.find(parameter)) != (parameters.end())) {
			return true;
		}else{
			cout << parameter << " is not a valid parameter in Mothur. Valid parameters are " << endl;
			for (it = parameters.begin(); it != parameters.end(); it++) {
				cout << it->first << ", ";
			}
			cout << endl;

			return false;
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidParameters class Function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidParameters class function isValidParameter. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

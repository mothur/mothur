/*
 *  optionparser.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "optionparser.h"

/***********************************************************************/

OptionParser::OptionParser(string option) {
	try {
		m = MothurOut::getInstance();
		if (option != "") {
			
			string key, value;		
			//reads in parameters and values
			while((option.find_first_of(',') != -1)) {  //while there are parameters
				m->splitAtComma(value, option);
				m->splitAtEquals(key, value);
				parameters[key] = value;
			}
			
			//in case there is no comma and to get last parameter after comma
			m->splitAtEquals(key, option);
			parameters[key] = option;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "OptionParser", "parse");
		exit(1);
	}
}

/***********************************************************************/

map<string, string> OptionParser::getParameters() {	return parameters;	}

/***********************************************************************/

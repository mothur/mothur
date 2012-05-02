#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

/*
 *  optionparser.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */



#include "mothur.h"
#include "mothurout.h"
#include "command.hpp"

/***********************************************************************/

class OptionParser {
public:
	OptionParser(string);
	~OptionParser() {}
	map<string, string> getParameters();
	bool getNameFile(vector<string>);
private:
	map<string, string> parameters;
	MothurOut* m;
};

/***********************************************************************/

#endif

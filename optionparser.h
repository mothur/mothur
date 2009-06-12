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


/***********************************************************************/

class OptionParser {
	public:
		OptionParser() {}
		~OptionParser() {}
		void parse(string, map<string, string>&);  //pass it an option string and a container
											  //fills the container key=parameter name, value=parameter value
};

/***********************************************************************/

#endif

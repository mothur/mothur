#ifndef CHIMERACOMMAND_H
#define CHIMERACOMMAND_H

/*
 *  chimeraseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"

/***********************************************************/

class ChimeraSeqsCommand : public Command {
public:
	ChimeraSeqsCommand(string);
	ChimeraSeqsCommand() {}
	~ChimeraSeqsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map< string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
		
private:
	vector<string> outputNames;
	map< string, vector<string> > outputTypes;

};

/***********************************************************/

#endif


#ifndef PHYLOTYPECOMMAND_H
#define PHYLOTYPECOMMAND_H


/*
 *  phylotypecommand.h
 *  Mothur
 *
 *  Created by westcott on 11/20/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "taxonomyequalizer.h"
#include "command.hpp"

/*************************************************************************/

class PhylotypeCommand : public Command {
	
public:
	PhylotypeCommand(string);	
	~PhylotypeCommand();
	int execute(); 
	void help();	
	
private:
	bool abort, allLines;
	string taxonomyFileName, label, outputDir;
	set<string> labels; //holds labels to be used
	int cutoff;
	
	map<int, int> currentNodes;
	map<int, int> parentNodes;
	map<int, int>::iterator itCurrent;

};


/*************************************************************************/


#endif




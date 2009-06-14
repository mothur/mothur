#ifndef READTREECOMMAND_H
#define READTREECOMMAND_H

/*
 *  readtreecommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "readtree.h"
#include "treemap.h"

class GlobalData;

class ReadTreeCommand : public Command {
public:
	ReadTreeCommand(string);
	~ReadTreeCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadTree* read;
	TreeMap* treeMap;
	string filename, treefile, groupfile;
	bool abort;

};


#endif

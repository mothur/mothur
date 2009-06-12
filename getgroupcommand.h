#ifndef GETGROUPCOMMAND_H
#define GETGROUPCOMMAND_H

/*
 *  getgroupcommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "globaldata.hpp"

class GetgroupCommand : public Command {
public:
	GetgroupCommand(string);
	~GetgroupCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	GroupMap* groupMap;
	string outputFile, sharedfile;
	ofstream out;
	ifstream in;
	bool abort;

};

#endif

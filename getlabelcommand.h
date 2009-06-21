#ifndef GETLABELCOMMAND_H
#define GETLABELCOMMAND_H

/*
 *  getlabelcommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "readmatrix.hpp"

class GlobalData;

class GetlabelCommand : public Command {
public:
	GetlabelCommand(string);
	~GetlabelCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	string filename;
	bool abort;
};

#endif

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
	GetlabelCommand();
	~GetlabelCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	string filename;
};

#endif
#ifndef GETGROUPCOMMAND_H
#define GETGROUPCOMMAND_H

/*
 *  getgroupcommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/2/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include "command.hpp"
#include "readmatrix.hpp"

class GlobalData;

class GetgroupCommand : public Command {
public:
	GetgroupCommand();
	~GetgroupCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	GroupMap* groupMap;
};

#endif

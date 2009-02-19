#ifndef GETLINECOMMAND_H
#define GETLINECOMMAND_H

/*
 *  getlinecommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "command.hpp"
#include "readmatrix.hpp"

class GlobalData;

class GetlineCommand : public Command {
public:
	GetlineCommand();
	~GetlineCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	string filename;
};

#endif

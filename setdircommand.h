#ifndef SETDIRCOMMAND_H
#define SETDIRCOMMAND_H

/*
 *  setoutdircommand.h
 *  Mothur
 *
 *  Created by westcott on 1/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "commandfactory.hpp"
#include "globaldata.hpp"

/**********************************************************/

class SetDirectoryCommand : public Command {
	
public:
	SetDirectoryCommand(string);
	~SetDirectoryCommand();
	int execute();
	void help();
	
private:
	CommandFactory* commandFactory;
	string output, input, tempdefault;
	bool abort;
		
};

/**********************************************************/
 
#endif


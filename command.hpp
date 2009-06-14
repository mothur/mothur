#ifndef COMMAND_HPP
#define COMMAND_HPP

/*
 *  command.h
 *  
 *
 *  Created by Pat Schloss on 10/23/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

/*This class is a parent to all the command classes.  */


#include "mothur.h"
#include "optionparser.h"
#include "validparameter.h"

class Command {
	public:
		virtual int execute() = 0;
		virtual void help() = 0;
};

#endif

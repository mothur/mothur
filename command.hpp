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
#include "mothurout.h"


class Command {
	
	public:
		Command() {  m = MothurOut::getInstance();   } 
		virtual vector<string> getValidParameters() = 0;
		virtual vector<string> getRequiredParameters() = 0; //adding "or" as the last element indicates one of the previous is needed
		virtual vector<string> getRequiredFiles() = 0; //adding "or" as the last element indicates one of the previous is needed
		virtual map<string, vector<string> > getOutputFiles() = 0; //file type to names
		virtual int execute() = 0;
		virtual void help() = 0;
		virtual ~Command() { }
	protected:
		MothurOut* m;
		bool calledHelp;
	
		map<string, vector<string> >::iterator itTypes;
};

#endif

#ifndef QUITCOMMAND_H
#define QUITCOMMAND_H
/*
 *  quitcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"

/* The quit() command:
	The quit command terminates the mothur program. 
	The quit command should be in the following format: quit ().   */


class QuitCommand : public Command {
	
public:
	QuitCommand(string);
	QuitCommand() {}
	~QuitCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	bool abort;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
};

#endif

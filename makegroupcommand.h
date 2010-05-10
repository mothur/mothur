#ifndef MAKEGROUPCOMMAND_H
#define MAKEGROUPCOMMAND_H

/*
 *  makegroupcommand.h
 *  Mothur
 *
 *  Created by westcott on 5/7/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"

class MakeGroupCommand : public Command {
	
public:
	MakeGroupCommand(string);	
	~MakeGroupCommand();
	int execute(); 
	void help();	
	
private:
		
	string fastaFileName, groups, outputDir, filename;
	vector<string> fastaFileNames;
	vector<string> groupsNames;
	
	bool abort;
};

#endif


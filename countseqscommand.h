#ifndef COuNTSEQSCOMMAND_H
#define COuNTSEQSCOMMAND_H

/*
 *  countseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 6/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"

class CountSeqsCommand : public Command {
	
public:
	
	CountSeqsCommand(string);
	CountSeqsCommand();	
	~CountSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "count.seqs";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Count.seqs"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	string namefile, groupfile, outputDir, groups;
	bool abort;
	vector<string> Groups;
};

#endif



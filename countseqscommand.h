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
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Count.seqs"; }
	string getDescription()		{ return "counts the number of sequences represented by each unique sequence in a namesfile"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	string namefile, groupfile, outputDir, groups;
	bool abort, large;
	vector<string> Groups, outputNames;
    int processors;
    
    int processSmall(string);
    int processLarge(string);
    map<int, string> processNameFile(string);
    map<int, string> getGroupNames(string, set<string>&);
    
};

#endif



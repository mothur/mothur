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
#include "groupmap.h"

class CountSeqsCommand : public Command {
	
public:
	
	CountSeqsCommand(string);
	~CountSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "count.seqs";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Count.seqs"; }
	string getDescription()		{ return "makes a count file from a names or shared file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
    
	string namefile, groupfile, outputDir, groups, sharedfile, countfile;
	bool abort, allLines, compress;
	vector<string> Groups, outputNames;
    set<string> labels;
    
    unsigned long long process(string);
    map<int, string> processNameFile(string);
    map<int, string> getGroupNames(string, set<string>&);
    
    unsigned long long processShared(vector<SharedRAbundVector*>& lookup, map<string, string> variables, vector<string>);
};

#endif



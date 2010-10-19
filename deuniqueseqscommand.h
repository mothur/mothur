#ifndef DEUNIQUESEQSCOMMAND_H
#define DEUNIQUESEQSCOMMAND_H
/*
 *  deuniqueseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"

/* This command is the reverse of unique.seqs */ 

class DeUniqueSeqsCommand : public Command {

public:
	DeUniqueSeqsCommand(string);
	DeUniqueSeqsCommand();
	~DeUniqueSeqsCommand() {}
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	
	
private:

	string fastaFile, nameFile, outputDir;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	bool abort;
	
	map<string, string> nameMap;
	
	int readNamesFile();
	
};

#endif


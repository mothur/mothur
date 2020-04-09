#ifndef MAKEFASTQCOMMAND_H
#define MAKEFASTQCOMMAND_H

/*
 *  makefastqcommand.h
 *  mothur
 *
 *  Created by westcott on 2/14/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"

class MakeFastQCommand : public Command {
	
public:
	
	MakeFastQCommand(string);	
	~MakeFastQCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "make.fastq";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Make.fastq"; }
	string getDescription()		{ return "creates a fastq file from a fasta and quality file"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	
	string fastafile, qualfile, outputDir, format;
	bool abort;
	vector<string> outputNames;
	
	string convertQual(vector<int>);
	
};

#endif



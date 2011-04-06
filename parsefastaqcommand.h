#ifndef PARSEFASTAQCOMMAND_H
#define PARSEFASTAQCOMMAND_H

/*
 *  parsefastaqcommand.h
 *  Mothur
 *
 *  Created by westcott on 9/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"

class ParseFastaQCommand : public Command {

public:
	ParseFastaQCommand(string);
	ParseFastaQCommand();
	~ParseFastaQCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "parse.fastq";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:

	vector<string> outputNames;	
	string outputDir, fastaQFile;
	bool abort;
	
	vector<int> convertQual(string);
};

#endif



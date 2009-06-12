#ifndef ALIGNCOMMAND_H
#define ALIGNCOMMAND_H

/*
 *  aligncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "globaldata.hpp"



class AlignCommand : public Command {
	
public:
	AlignCommand(string);	
	~AlignCommand();
	int execute(); 
	void help();	

private:
	GlobalData* globaldata;
	OptionParser* parser;
	map<string, string> parameters;
	map<string, string>::iterator it;
	bool abort;
	string candidateFileName, templateFileName, distanceFileName, search, align;
	int kmerSize;
	float match, misMatch, gapOpen, gapExtend;
	ofstream out;
	ifstream in;
	int ableToOpen;
	

};



#endif
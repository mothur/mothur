#ifndef FILTERSEQSCOMMAND_H
#define FILTERSEQSCOMMAND_H

/*
 *  filterseqscommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "globaldata.hpp"
#include "sequence.hpp"

class FilterSeqsCommand : public Command {

public:
	FilterSeqsCommand(string);
	~FilterSeqsCommand() {};
	int execute();	
	void help();
	
private:
	void doHard();
	void doTrump(Sequence);
	void doVertical();
	void doSoft();
	void getFreqs(Sequence);
	string vertical, filter, fastafile, hard;	
	int alignmentLength;

	char trump;
	bool abort;
	float soft;
	int numSeqs;
	OptionParser* parser;
	map<string, string> parameters;
	map<string, string>::iterator it;
	
	GlobalData* globaldata;	
	vector<int> a, t, g, c, gap;

};

#endif

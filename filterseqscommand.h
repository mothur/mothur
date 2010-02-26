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
#include "filters.h"

class Sequence;
class FilterSeqsCommand : public Command {

public:
	FilterSeqsCommand(string);
	~FilterSeqsCommand() {};
	int execute();	
	void help();
	
private:
	string vertical, filter, fasta, hard, outputDir, filterFileName;
	vector<string> fastafileNames;	
	int alignmentLength;

	char trump;
	bool abort;
	float soft;
	int numSeqs;
	
	Filters F;
		
	vector<int> a, t, g, c, gap;

};

#endif

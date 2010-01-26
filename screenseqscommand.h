#ifndef SCREENSEQSCOMMAND_H
#define SCREENSEQSCOMMAND_H

/*
 *  screenseqscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/3/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */
#include "mothur.h"
#include "command.hpp"

class ScreenSeqsCommand : public Command {
	
public:
	ScreenSeqsCommand(string);
	~ScreenSeqsCommand();
	int execute();
	void help();
	
private:
	void screenNameGroupFile(set<string>);
	void screenGroupFile(set<string>);
	void screenAlignReport(set<string>);
	
	bool abort;
	string fastafile, namefile, groupfile, alignreport, outputDir;
	int startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength;
};

#endif

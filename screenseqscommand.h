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
#include "globaldata.hpp"
#include "readfasta.h"
#include "readnexus.h"
#include "readclustal.h"
#include "readseqsphylip.h"

using namespace std;

class ScreenSeqsCommand : public Command {
	
public:
	ScreenSeqsCommand();
	~ScreenSeqsCommand();
	int execute();
private:
	void screenNameGroupFile(set<string>);
	void screenGroupFile(set<string>);

	GlobalData* globaldata;	
};

#endif

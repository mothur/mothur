#ifndef REVERSECOMMAND_H
#define REVERSECOMMAND_H

/*
 *  reversecommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"

class ReverseSeqsCommand : public Command {
public:
	ReverseSeqsCommand(string);
	~ReverseSeqsCommand();
	int execute();
	void help();
	
private:

	bool abort;
	string fasta;
	
};

#endif

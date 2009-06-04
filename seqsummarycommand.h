#ifndef SEQSUMMARYCOMMAND_H
#define SEQSUMMARYCOMMAND_H

/*
 *  seqcoordcommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "globaldata.hpp"



class SeqSummaryCommand : public Command {
public:
	SeqSummaryCommand();
	~SeqSummaryCommand();
	int execute();
	
private:
	GlobalData* globaldata;	

};

#endif

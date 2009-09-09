#ifndef CHIMERACOMMAND_H
#define CHIMERACOMMAND_H

/*
 *  chimeraseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "chimera.h"


/***********************************************************/

class ChimeraSeqsCommand : public Command {
public:
	ChimeraSeqsCommand(string);
	~ChimeraSeqsCommand();
	int execute();
	void help();
	
		
private:
	
	bool abort;
	string method, fastafile, templatefile, consfile, quanfile, maskfile;
	bool filter, correction;
	int processors, midpoint, averageLeft, averageRight, window, iters, increment, numwanted, ksize;
	Chimera* chimera;
	
	
};

/***********************************************************/

#endif


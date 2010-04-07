/*
 *  chimeraseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimeraseqscommand.h"

//***************************************************************************************************************

ChimeraSeqsCommand::ChimeraSeqsCommand(string option)  {}
//**********************************************************************************************************************

void ChimeraSeqsCommand::help(){}

//***************************************************************************************************************

ChimeraSeqsCommand::~ChimeraSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraSeqsCommand::execute(){
	
		m->mothurOut("The chimera.seqs command has been broken up into 5 separate commands.\n");
		m->mothurOut("The chimera.bellerophon, chimera.ccode, chimera.check, chimera.pintail and chimera.slayer commands.\n");
	
	return 0;
}
//**********************************************************************************************************************



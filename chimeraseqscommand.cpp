/*
 *  chimeraseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimeraseqscommand.h"

//**********************************************************************************************************************
vector<string> ChimeraSeqsCommand::getValidParameters(){	
	try {
		vector<string> myArray; 
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ChimeraSeqsCommand::getRequiredParameters(){	
	try {
		vector<string> myArray; 
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ChimeraSeqsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "getRequiredFiles");
		exit(1);
	}
}
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



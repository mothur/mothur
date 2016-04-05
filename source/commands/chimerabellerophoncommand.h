#ifndef CHIMERABELLEROPHONCOMMAND_H
#define CHIMERABELLEROPHONCOMMAND_H

/*
 *  chimerabellerophoncommand.h
 *  Mothur
 *
 *  Created by westcott on 4/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "mothurchimera.h"


/***********************************************************/

class ChimeraBellerophonCommand : public Command {
public:
	ChimeraBellerophonCommand(string);
	ChimeraBellerophonCommand();
	~ChimeraBellerophonCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.bellerophon";	}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Huber T, Faulkner G, Hugenholtz P (2004). Bellerophon: a program to detect chimeric sequences in multiple sequence alignments. Bioinformatics 20: 2317-9. \nhttp://www.mothur.org/wiki/Chimera.bellerophon"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
		
private:

	bool abort, filter, correction;
	string fastafile, outputDir;
	int processors, window, increment, numSeqs;
	MothurChimera* chimera;
	vector<string> outputNames;
	vector<string> fastaFileNames;
};

/***********************************************************/

#endif



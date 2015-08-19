#ifndef CHIMERACHECKCOMMAND_H
#define CHIMERACHECKCOMMAND_H

/*
 *  chimeracheckcommand.h
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "chimera.h"
#include "chimeracheckrdp.h"


/***********************************************************/

class ChimeraCheckCommand : public Command {
public:
	ChimeraCheckCommand(string);
	ChimeraCheckCommand();
	~ChimeraCheckCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.check";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "CHIMERA_CHECK version 2.7 written by Niels Larsen (http://wdcm.nig.ac.jp/RDP/docs/chimera_doc.html) \nhttp://www.mothur.org/wiki/Chimera.check"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	
	int driver(linePair*, string, string);
	int createProcesses(string, string);

	bool abort, svg, save;
	string fastafile, templatefile, namefile, outputDir;
	int processors, increment, ksize, numSeqs, templateSeqsLength;
	Chimera* chimera;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	vector<string> outputNames;
};

/***********************************************************/

#endif



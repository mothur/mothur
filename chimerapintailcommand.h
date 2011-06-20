#ifndef CHIMERAPINTAILCOMMAND_H
#define CHIMERAPINTAILCOMMAND_H

/*
 *  chimerapintailcommand.h
 *  Mothur
 *
 *  Created by westcott on 4/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "chimera.h"


/***********************************************************/

class ChimeraPintailCommand : public Command {

public:

	ChimeraPintailCommand(string);
	ChimeraPintailCommand();
	~ChimeraPintailCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.pintail";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "Ashelford KE, Chuzhanova NA, Fry JC, Jones AJ, Weightman AJ (2005). At least 1 in 20 16S rRNA sequence records currently held in public repositories is estimated to contain substantial anomalies. Appl Environ Microbiol 71: 7724-36. \nAshelford KE, Chuzhanova NA, Fry JC, Jones AJ, Weightman AJ (2006). New screening software shows that most recent large 16S rRNA gene clone libraries contain chimeras. Appl Environ Microbiol 72: 5734-41. \nhttp://www.mothur.org/wiki/Chimera.pintail"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:

	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};

	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	
	int driver(linePair*, string, string, string);
	int createProcesses(string, string, string);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long int>&);
	#endif

	bool abort, filter;
	string fastafile, templatefile, consfile, quanfile, maskfile, outputDir, inputDir;
	int processors, window, increment, numSeqs, templateSeqsLength;
	Chimera* chimera;
	vector<string> outputNames;
	vector<string> fastaFileNames;
};

/***********************************************************/

#endif



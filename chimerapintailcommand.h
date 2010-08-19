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
	~ChimeraPintailCommand();
	int execute();
	void help();
	
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



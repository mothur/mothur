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
		int start;
		int numSeqs;
		linePair(long int i, int j) : start(i), numSeqs(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	
	int driver(linePair*, string, string, string);
	int createProcesses(string, string, string);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, vector<long>&);
	#endif

	bool abort, filter, MPIWroteAccnos;
	string fastafile, templatefile, consfile, quanfile, maskfile, outputDir, inputDir;
	int processors, window, increment, numSeqs, templateSeqsLength;
	Chimera* chimera;
	vector<string> outputNames;
	vector<string> fastaFileNames;
	
	
};

/***********************************************************/

#endif



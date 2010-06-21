#ifndef CHIMERASLAYERCOMMAND_H
#define CHIMERASLAYERCOMMAND_H

/*
 *  chimeraslayercommand.h
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "chimera.h"


/***********************************************************/

class ChimeraSlayerCommand : public Command {
public:
	ChimeraSlayerCommand(string);
	~ChimeraSlayerCommand();
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

	bool abort, realign, MPIWroteAccnos;
	string fastafile, templatefile, outputDir, search;
	int processors, window, iters, increment, numwanted, ksize, match, mismatch, parents, minSimilarity, minCoverage, minBS, minSNP, numSeqs, templateSeqsLength;
	float divR;
	Chimera* chimera;
	
	vector<string> outputNames;
	vector<string> fastaFileNames;
	
};

/***********************************************************/

#endif



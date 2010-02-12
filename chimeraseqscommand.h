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

	struct linePair {
		int start;
		int numSeqs;
		linePair(long int i, int j) : start(i), numSeqs(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	
	int driver(linePair*, string, string);
	void createProcesses(string, string);
	void appendOutputFiles(string, string); 

	bool abort;
	string method, fastafile, templatefile, consfile, quanfile, maskfile, namefile, outputDir;
	bool filter, correction, svg, printAll;
	int processors, midpoint, averageLeft, averageRight, window, iters, increment, numwanted, ksize, match, mismatch, parents, minSimilarity, minCoverage, minBS, minSNP, numSeqs;
	float divR;
	Chimera* chimera;
	
	
};

/***********************************************************/

#endif

